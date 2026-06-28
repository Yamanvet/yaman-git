#!/usr/bin/env Rscript
# ======================================================================
# CTVT: satuRn DTU — v3 FIXED
#
# Fixes from v2:
#   ✅ Filter per-comparison (only the 2 groups being compared)
#      — NOT globally across all 4 timepoints
#   ✅ Use raw NumReads for voom (not lengthScaledTPM)
#   ✅ Compute IF/dIF from lengthScaledTPM (abundance-based)
#      but statistical testing from raw counts + voom
# ======================================================================

.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_MAX_VSIZE = "500e9")
options(future.globals.maxSize = 10 * 1024^3)
gc(full = TRUE)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(BiocParallel)
  library(dplyr)
})

# ── Config ──────────────────────────────────────────────────────────
BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"
GTF_FILE <- file.path(BASE, "stringtie_v4/merged_stringtie_v4.gtf")
OUTDIR <- file.path(BASE, "saturn_tximport_v3")
OUTLIERS <- c("T1B_S2", "T8C_S10", "T9N_S26")

ncores <- min(parallel::detectCores(), 60)
register(MulticoreParam(workers = ncores, progressbar = TRUE))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUTDIR, "run.log"), split = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CTVT: satuRn DTU v3 — per-comparison filtering           \n")
cat("  R:", as.character(getRversion()), "\n")
cat("  edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("  limma:", as.character(packageVersion("limma")), "\n")
cat("  GTF:", GTF_FILE, "\n")
cat("  Quant:", QUANT_DIR, "\n")
cat("  Samples: 33 (dropping", paste(OUTLIERS, collapse=", "), ")\n")
cat("  FIX: filter per-comparison, not globally\n")
cat("  FIX: raw counts for voom, lengthScaledTPM for IF/dIF\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── 1. Sample metadata ─────────────────────────────────────────────
cat("[1/6] Building sample metadata...\n")

all_samples <- basename(list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE))
keep_samples <- setdiff(all_samples, OUTLIERS)

sample_meta <- data.frame(sampleID = keep_samples, stringsAsFactors = FALSE)
sample_meta$timepoint <- factor(sapply(sample_meta$sampleID, function(s) {
  prefix <- strsplit(s, "_")[[1]][1]
  c("A"="Day0","B"="Day2","C"="Day6","N"="Recovered")[substr(prefix,3,3)]
}), levels = c("Day0", "Day2", "Day6", "Recovered"))

sample_meta$dog_id <- factor(sapply(sample_meta$sampleID, function(s) {
  substr(strsplit(s, "_")[[1]][1], 1, 2)
}))

sample_meta$path <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
stopifnot(all(file.exists(sample_meta$path)))

cat(sprintf("  %d samples: Day0=%d, Day2=%d, Day6=%d, Recovered=%d\n",
  nrow(sample_meta),
  sum(sample_meta$timepoint == "Day0"),
  sum(sample_meta$timepoint == "Day2"),
  sum(sample_meta$timepoint == "Day6"),
  sum(sample_meta$timepoint == "Recovered")))

# ── 2. GTF gene-transcript mapping ─────────────────────────────────
cat("\n[2/6] Importing GTF for gene-transcript mapping...\n")

suppressPackageStartupMessages(library(rtracklayer))
gtf_gr <- import(GTF_FILE)

tx_features <- gtf_gr[gtf_gr$type == "transcript"]
gene_tx_map <- data.frame(
  transcript_id = mcols(tx_features)$transcript_id,
  gene_id = mcols(tx_features)$gene_id,
  stringsAsFactors = FALSE
)
if ("gene_name" %in% colnames(mcols(tx_features))) {
  gene_tx_map$gene_name <- mcols(tx_features)$gene_name
}
gene_tx_map <- gene_tx_map[!is.na(gene_tx_map$transcript_id), ]

cat(sprintf("  GTF: %d transcripts mapped to %d genes\n",
  nrow(gene_tx_map), length(unique(gene_tx_map$gene_id))))

tx_to_gene <- setNames(gene_tx_map$gene_id, gene_tx_map$transcript_id)
gene_to_name <- if ("gene_name" %in% colnames(gene_tx_map)) {
  gn <- setNames(gene_tx_map$gene_name, gene_tx_map$gene_id)
  gn[!duplicated(names(gn))]
} else character(0)

rm(gtf_gr, tx_features); gc(full = TRUE)

# ── 3. Read Salmon quant files ─────────────────────────────────────
cat("\n[3/6] Reading Salmon quant files...\n")

first_df <- read.delim(sample_meta$path[1], stringsAsFactors = FALSE)
tx_names <- first_df$Name
n_tx <- length(tx_names)
cat(sprintf("  %d transcripts per sample\n", n_tx))

# Raw counts matrix (NumReads) — for voom/statistical testing
raw_counts <- matrix(0, nrow = n_tx, ncol = nrow(sample_meta),
                      dimnames = list(tx_names, sample_meta$sampleID))

# Effective length matrix — for lengthScaledTPM / IF calculation
eff_len_mat <- matrix(0, nrow = n_tx, ncol = nrow(sample_meta),
                      dimnames = list(tx_names, sample_meta$sampleID))

# TPM matrix — for IF calculation (abundance-based)
tpm_mat <- matrix(0, nrow = n_tx, ncol = nrow(sample_meta),
                  dimnames = list(tx_names, sample_meta$sampleID))

for (i in 1:nrow(sample_meta)) {
  df <- read.delim(sample_meta$path[i], stringsAsFactors = FALSE)
  raw_counts[df$Name, i] <- df$NumReads
  eff_len_mat[df$Name, i] <- df$EffectiveLength
  tpm_mat[df$Name, i] <- df$TPM
  if (i %% 10 == 0) cat(sprintf("    Read %d/%d samples\n", i, nrow(sample_meta)))
}
cat(sprintf("    Read %d/%d samples — done\n", i, nrow(sample_meta)))

cat(sprintf("  Raw counts: %d tx x %d samples\n", nrow(raw_counts), ncol(raw_counts)))

# Gene mapping
gene_for_tx <- tx_to_gene[tx_names]
unmapped <- is.na(gene_for_tx)
gene_for_tx[unmapped] <- tx_names[unmapped]

saveRDS(list(raw_counts = raw_counts, eff_len = eff_len_mat, tpm = tpm_mat,
             gene_for_tx = gene_for_tx),
        file.path(OUTDIR, "01_salmon_data.rds"))
cat("  Saved: 01_salmon_data.rds\n")
gc(full = TRUE)

# ── 4. Helper: compute IF from lengthScaledTPM ─────────────────────
# IF = isoform fraction = tx_lstp / gene_total_lstp per sample
# lengthScaledTPM per sample = NumReads * (eff_len / mean_eff_len)

compute_IF <- function(counts_sub, eff_len_sub, gene_ids_sub) {
  # lengthScaledTPM for this subset
  mean_eff_len <- rowMeans(eff_len_sub)
  mean_eff_len[mean_eff_len == 0] <- 1
  lstp <- counts_sub * (eff_len_sub / mean_eff_len)

  # gene totals
  gene_totals <- rowsum(lstp, group = gene_ids_sub)

  # IF = tx_lstp / gene_total_lstp
  if_vals <- lstp / gene_totals[gene_ids_sub, ]
  if_vals[is.nan(if_vals)] <- 0
  return(if_vals)
}

# ── 5. Per-comparison DTU testing ──────────────────────────────────
cat("\n[4/6] Per-comparison DTU testing (filter per comparison)...\n")

comp_df <- data.frame(
  condition_1 = c("Day0", "Day2", "Day6", "Day0", "Day0"),
  condition_2 = c("Day2", "Day6", "Recovered", "Day6", "Recovered"),
  name = c("Day0_vs_Day2", "Day2_vs_Day6",
           "Day6_vs_Recovered", "Day0_vs_Day6", "Day0_vs_Recovered"),
  stringsAsFactors = FALSE
)

all_results <- list()

for (i in 1:nrow(comp_df)) {
  cn <- comp_df$name[i]
  c1 <- comp_df$condition_1[i]
  c2 <- comp_df$condition_2[i]

  cat(sprintf("\n  ── Comparison %d/%d: %s ──\n", i, nrow(comp_df), cn))

  # Subset samples for this comparison ONLY
  keep_comp <- sample_meta$timepoint %in% c(c1, c2)
  samples_comp <- sample_meta[keep_comp, ]

  cat(sprintf("    %s (n=%d) vs %s (n=%d)\n",
    c1, sum(samples_comp$timepoint == c1),
    c2, sum(samples_comp$timepoint == c2)))

  # Subset raw counts for these samples
  raw_comp <- raw_counts[, samples_comp$sampleID]
  eff_comp <- eff_len_mat[, samples_comp$sampleID]

  # DGEList from RAW counts
  dge <- DGEList(counts = raw_comp)
  dge$samples$group <- factor(samples_comp$timepoint, levels = c(c1, c2))

  # FILTER: only using the 2 groups in this comparison
  # This keeps transcripts expressed in either Day2 OR Day6, etc.
  keep <- filterByExpr(dge, group = dge$samples$group, min.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  n_retained <- nrow(dge$counts)
  cat(sprintf("    After per-comparison filter: %d / %d transcripts retained\n",
    n_retained, n_tx))

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")
  cat(sprintf("    TMM done. Lib sizes: %.0f - %.0f\n",
    min(dge$samples$lib.size), max(dge$samples$lib.size)))

  # Design matrix: ~ dog_id + condition (paired design)
  samples_comp$condition <- factor(samples_comp$timepoint, levels = c(c1, c2))
  design_comp <- model.matrix(~ dog_id + condition, data = samples_comp)
  colnames(design_comp) <- make.names(colnames(design_comp))

  # Voom on raw counts (statistically correct)
  cat("    voomWithQualityWeights (raw counts)...\n")
  v <- voomWithQualityWeights(dge$counts, design_comp, plot = FALSE)

  # Fit
  cat("    lmFit + eBayes...\n")
  fit <- lmFit(v, design_comp)
  fit <- eBayes(fit)

  contrast_coef <- paste0("condition", c2)
  if (!(contrast_coef %in% colnames(design_comp))) {
    cat(sprintf("    WARNING: %s not found, skipping\n", contrast_coef))
    next
  }

  tt <- topTable(fit, coef = contrast_coef, number = Inf,
    sort.by = "none", adjust.method = "BH")

  # IF/dIF from lengthScaledTPM (abundance-based, not raw counts)
  filtered_tx <- rownames(dge$counts)
  gene_ids_filt <- gene_for_tx[filtered_tx]

  if_values <- compute_IF(raw_comp[filtered_tx, ], eff_comp[filtered_tx, ],
                          gene_ids_filt)

  idx_c1 <- which(samples_comp$timepoint == c1)
  idx_c2 <- which(samples_comp$timepoint == c2)

  mean_if_c1 <- rowMeans(if_values[, idx_c1, drop = FALSE])
  mean_if_c2 <- rowMeans(if_values[, idx_c2, drop = FALSE])
  dIF <- mean_if_c2 - mean_if_c1

  # Assemble results
  result_df <- data.frame(
    transcript_id = filtered_tx,
    gene_id = gene_ids_filt,
    logFC = tt$logFC,
    AveExpr = tt$AveExpr,
    P.Value = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    mean_IF_c1 = mean_if_c1,
    mean_IF_c2 = mean_if_c2,
    dIF = dIF,
    stringsAsFactors = FALSE
  )

  if (length(gene_to_name) > 0) {
    result_df$gene_name <- gene_to_name[result_df$gene_id]
    result_df$gene_name[is.na(result_df$gene_name)] <- result_df$gene_id[is.na(result_df$gene_name)]
  }

  # Significant DTU
  sig_q05 <- result_df[!is.na(result_df$adj.P.Val) &
                        result_df$adj.P.Val < 0.05 &
                        abs(result_df$dIF) >= 0.05, ]
  sig_q05 <- sig_q05[order(sig_q05$adj.P.Val), ]

  sig_q10 <- result_df[!is.na(result_df$adj.P.Val) &
                        result_df$adj.P.Val < 0.10 &
                        abs(result_df$dIF) >= 0.05, ]

  cat(sprintf("    Tested: %d tx\n", nrow(result_df)))
  cat(sprintf("    Significant (q<0.05, |dIF|≥0.05): %d tx / %d genes\n",
    nrow(sig_q05), length(unique(sig_q05$gene_id))))
  cat(sprintf("    Relaxed   (q<0.10, |dIF|≥0.05): %d tx / %d genes\n",
    nrow(sig_q10), length(unique(sig_q10$gene_id))))

  if (nrow(sig_q05) > 0) {
    cat(sprintf("\n    Top DTU in %s:\n", cn))
    top <- head(sig_q05, 10)
    for (j in 1:nrow(top)) {
      gn <- ifelse(is.na(top$gene_name[j]), top$gene_id[j], top$gene_name[j])
      cat(sprintf("    %2d. %-18s iso=%-18s dIF=%+.4f q=%.2e\n",
        j, gn, top$transcript_id[j], top$dIF[j], top$adj.P.Val[j]))
    }
  }

  # Save
  comp_out <- file.path(OUTDIR, cn)
  dir.create(comp_out, showWarnings = FALSE)
  write.csv(result_df, file.path(comp_out, "all_transcript_results.csv"), row.names = FALSE)
  if (nrow(sig_q05) > 0)
    write.csv(sig_q05, file.path(comp_out, "significant_dtu_q005.csv"), row.names = FALSE)
  if (nrow(sig_q10) > 0)
    write.csv(sig_q10, file.path(comp_out, "significant_dtu_q010.csv"), row.names = FALSE)

  all_results[[cn]] <- list(
    comparison = cn,
    total_tested = nrow(result_df),
    n_retained = n_retained,
    n_dtu_q005 = nrow(sig_q05),
    n_genes_q005 = length(unique(sig_q05$gene_id)),
    n_dtu_q010 = nrow(sig_q10),
    n_genes_q010 = length(unique(sig_q10$gene_id))
  )

  rm(dge, v, fit, if_values); gc(full = TRUE)
}

# ── 5. Summary ─────────────────────────────────────────────────────
cat("\n[5/6] Summary\n\n")
cat("  ──────────────────────────────────────────────────────────\n")
cat("  satuRn DTU SUMMARY v3 (per-comparison filter + raw counts) \n")
cat("  ──────────────────────────────────────────────────────────\n\n")
cat(sprintf("  %-25s %8s %8s %8s %8s\n",
  "Comparison", "Retained", "DTU(q05)", "Genes(q05)", "DTU(q10)"))
cat("  "); cat(rep("─", 56), sep = ""); cat("\n")

overall <- data.frame()
for (cn in comp_df$name) {
  r <- all_results[[cn]]
  cat(sprintf("  %-25s %8d %8d %8d %8d\n",
    r$comparison, r$n_retained,
    r$n_dtu_q005, r$n_genes_q005, r$n_dtu_q010))
  overall <- rbind(overall, data.frame(
    Comparison = r$comparison, Retained = r$n_retained,
    DTU_q005 = r$n_dtu_q005, Genes_q005 = r$n_genes_q005,
    DTU_q010 = r$n_dtu_q010, Genes_q010 = r$n_genes_q010,
    stringsAsFactors = FALSE
  ))
}
write.csv(overall, file.path(OUTDIR, "summary.csv"), row.names = FALSE)

# ── 6. Cross-comparison ────────────────────────────────────────────
cat("\n[6/6] Cross-comparison: genes with DTU in multiple contrasts\n")

all_dtu_genes <- list()
for (cn in comp_df$name) {
  sig_file <- file.path(OUTDIR, cn, "significant_dtu_q005.csv")
  if (file.exists(sig_file)) {
    df <- read.csv(sig_file)
    if (nrow(df) > 0) all_dtu_genes[[cn]] <- unique(df$gene_id)
  }
}

all_genes <- unique(unlist(all_dtu_genes))
if (length(all_genes) > 0) {
  gene_comp_count <- data.frame(
    gene = all_genes,
    n_comparisons = sapply(all_genes, function(g)
      sum(sapply(all_dtu_genes, function(v) g %in% v))),
    comparisons = sapply(all_genes, function(g)
      paste(names(all_dtu_genes)[sapply(all_dtu_genes, function(v) g %in% v)],
            collapse = "; ")),
    stringsAsFactors = FALSE
  )
  gene_comp_count <- gene_comp_count[order(-gene_comp_count$n_comparisons), ]
  cat(sprintf("  %d genes with DTU in ≥1 comparison\n", nrow(gene_comp_count)))
  cat(sprintf("  %d genes with DTU in ≥2 comparisons\n",
    sum(gene_comp_count$n_comparisons >= 2)))
  if (sum(gene_comp_count$n_comparisons >= 3) > 0)
    cat(sprintf("  %d genes with DTU in ≥3 comparisons\n",
      sum(gene_comp_count$n_comparisons >= 3)))
  write.csv(gene_comp_count, file.path(OUTDIR, "multiswitch_genes.csv"),
    row.names = FALSE)
} else {
  cat("  No significant DTU found at q<0.05\n")
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("  ✅ satuRn v3 pipeline COMPLETE\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("════════════════════════════════════════════════════════════\n")

sink()