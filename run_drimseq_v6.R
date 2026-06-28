#!/usr/bin/env Rscript
# ======================================================================
# CTVT: DRIMSeq DTU v6 — Beta-binomial GLM (less conservative than satuRn)
#
# DRIMSeq uses beta-binomial GLM on transcript proportions within genes.
# Standard B-H FDR correction (not empirical null like satuRn).
#
# Fixes from v6 initial attempt:
#   - Use dmPrecision/dmFit/dmTest (current DRIMSeq API)
#   - geneIDs() is not exported; use counts(d) to get gene_id column
#   - Proper design matrix with model.matrix
#   - dmFilter with correct parameters for small n
# ======================================================================

.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_MAX_VSIZE = "500e9")
options(future.globals.maxSize = 10 * 1024^3)
gc(full = TRUE)

suppressPackageStartupMessages({
  library(DRIMSeq)
  library(edgeR)
})

# ── Config ──────────────────────────────────────────────────────────
BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"
GTF_FILE <- file.path(BASE, "stringtie_v4/merged_stringtie_v4.gtf")
OUTDIR <- file.path(BASE, "drimseq_v6")
OUTLIERS <- c("T1B_S2", "T8C_S10", "T9N_S26")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUTDIR, "run.log"), split = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CTVT: DRIMSeq DTU v6 — beta-binomial GLM\n")
cat("  R:", as.character(getRversion()), "\n")
cat("  DRIMSeq:", as.character(packageVersion("DRIMSeq")), "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── 1. Sample metadata ─────────────────────────────────────────────
cat("[1/5] Building sample metadata...\n")

all_samples <- basename(list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE))
keep_samples <- setdiff(all_samples, OUTLIERS)

sample_meta <- data.frame(sampleID = keep_samples, stringsAsFactors = FALSE)
sample_meta$timepoint <- factor(sapply(sample_meta$sampleID, function(s) {
  prefix <- strsplit(s, "_")[[1]][1]
  c("A"="Day0","B"="Day2","C"="Day6","N"="Recovered")[substr(prefix,3,3)]
}), levels = c("Day0", "Day2", "Day6", "Recovered"))

sample_meta$path <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
stopifnot(all(file.exists(sample_meta$path)))

cat(sprintf("  %d samples: Day0=%d, Day2=%d, Day6=%d, Recovered=%d\n",
  nrow(sample_meta),
  sum(sample_meta$timepoint == "Day0"),
  sum(sample_meta$timepoint == "Day2"),
  sum(sample_meta$timepoint == "Day6"),
  sum(sample_meta$timepoint == "Recovered")))

# ── 2. GTF gene-transcript mapping ─────────────────────────────────
cat("\n[2/5] Importing GTF for gene-transcript mapping...\n")

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

# ── 3. Load Salmon counts ───────────────────────────────────────────
cat("\n[3/5] Loading Salmon counts...\n")

cached <- file.path("/home/vet/CTVT_raw_fastq_36samples/saturn_v4/01_salmon_counts.rds")
if (file.exists(cached)) {
  loaded <- readRDS(cached)
  raw_counts <- loaded$raw_counts
  gene_for_tx <- loaded$gene_for_tx
  cat(sprintf("  Loaded %d tx x %d samples from cache\n",
    nrow(raw_counts), ncol(raw_counts)))
} else {
  first_df <- read.delim(sample_meta$path[1], stringsAsFactors = FALSE)
  tx_names <- first_df$Name
  n_tx <- length(tx_names)
  raw_counts <- matrix(0, nrow = n_tx, ncol = nrow(sample_meta),
                        dimnames = list(tx_names, sample_meta$sampleID))
  for (i in 1:nrow(sample_meta)) {
    df <- read.delim(sample_meta$path[i], stringsAsFactors = FALSE)
    raw_counts[df$Name, i] <- df$NumReads
  }
  gene_for_tx <- tx_to_gene[tx_names]
  unmapped <- is.na(gene_for_tx)
  gene_for_tx[unmapped] <- tx_names[unmapped]
}
gc(full = TRUE)

# ── 4. Per-comparison DRIMSeq DTU ──────────────────────────────────
cat("\n[4/5] Per-comparison DRIMSeq DTU (beta-binomial GLM)...\n")

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

  cat(sprintf("\n  -- Comparison %d/%d: %s --\n", i, nrow(comp_df), cn))

  keep_comp <- sample_meta$timepoint %in% c(c1, c2)
  samples_comp <- sample_meta[keep_comp, ]

  cat(sprintf("    %s (n=%d) vs %s (n=%d)\n",
    c1, sum(samples_comp$timepoint == c1),
    c2, sum(samples_comp$timepoint == c2)))

  raw_comp <- raw_counts[, samples_comp$sampleID]

  # DGEList for filtering
  dge <- DGEList(counts = raw_comp)
  dge$samples$group <- factor(samples_comp$timepoint, levels = c(c1, c2))

  keep <- filterByExpr(dge, group = dge$samples$group, min.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  filtered_tx <- rownames(dge$counts)
  gene_ids_filt <- gene_for_tx[filtered_tx]

  n_retained <- length(filtered_tx)
  cat(sprintf("    After per-comparison filter: %d tx\n", n_retained))

  # Multi-isoform genes only
  gene_counts <- table(gene_ids_filt)
  multi_genes <- names(gene_counts[gene_counts >= 2])
  keep_multi <- gene_ids_filt %in% multi_genes

  n_multi <- sum(keep_multi)
  n_multi_genes <- length(multi_genes)
  cat(sprintf("    Multi-isoform: %d genes (%d tx)\n", n_multi_genes, n_multi))

  if (n_multi == 0) { next }

  dge_multi <- dge[keep_multi, ]
  tx_multi <- rownames(dge_multi$counts)
  gene_multi <- gene_for_tx[tx_multi]

  # ── Build DRIMSeq data frame (feature_id, gene_id, sample counts) ──
  counts_df <- as.data.frame(dge_multi$counts)
  counts_df$feature_id <- tx_multi
  counts_df$gene_id <- gene_multi

  # Reorder columns: feature_id, gene_id, then samples
  counts_df <- counts_df[, c("feature_id", "gene_id", samples_comp$sampleID)]

  sample_info <- data.frame(
    sample_id = samples_comp$sampleID,
    condition = factor(samples_comp$timepoint, levels = c(c1, c2)),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    DRIMSeq input: %d tx x %d samples\n", n_multi, nrow(sample_info)))

  # Create dmDSdata object
  d <- dmDSdata(counts = counts_df, samples = sample_info)

  cat(sprintf("    dmDSdata created: %d genes, %d features\n",
    length(unique(counts(d)$gene_id)), nrow(counts(d))))

  # Filter low-expression transcripts
  n_min_samps <- min(sum(samples_comp$timepoint == c1),
                     sum(samples_comp$timepoint == c2))
  cat(sprintf("    dmFilter (min_samps=%d)...\n", n_min_samps))
  d <- dmFilter(d,
    min_samps_gene_expr = n_min_samps,
    min_samps_feature_expr = n_min_samps,
    min_samps_feature_prop = n_min_samps,
    min_gene_expr = 10,
    min_feature_expr = 5,
    min_feature_prop = 0.01)

  n_genes_after <- length(unique(counts(d)$gene_id))
  n_tx_after <- nrow(counts(d))
  cat(sprintf("    After dmFilter: %d genes, %d tx\n", n_genes_after, n_tx_after))

  # Design matrix
  design_full <- model.matrix(~ condition, data = samples(d))

  # Estimate precision (beta-binomial dispersion)
  cat("    dmPrecision...\n")
  set.seed(42)
  d <- dmPrecision(d, design = design_full)

  cat(sprintf("    Common precision: %.4f\n", common_precision(d)))

  # Fit model
  cat("    dmFit...\n")
  d <- dmFit(d, design = design_full, verbose = 1)

  # Test for DTU
  cat("    dmTest...\n")
  d <- dmTest(d, coef = "condition", verbose = 1)

  # Extract results
  res_gene <- results(d)
  res_tx <- results(d, level = "feature")

  cat(sprintf("    Genes tested: %d\n", nrow(res_gene)))
  cat(sprintf("    Transcripts tested: %d\n", nrow(res_tx)))

  # Gene-level DTU
  sig_genes_q05 <- res_gene[!is.na(res_gene$adj_pvalue) & res_gene$adj_pvalue < 0.05, ]
  sig_genes_q10 <- res_gene[!is.na(res_gene$adj_pvalue) & res_gene$adj_pvalue < 0.10, ]
  sig_genes_q20 <- res_gene[!is.na(res_gene$adj_pvalue) & res_gene$adj_pvalue < 0.20, ]

  cat(sprintf("    Gene-level: q<0.05=%d, q<0.10=%d, q<0.20=%d\n",
    nrow(sig_genes_q05), nrow(sig_genes_q10), nrow(sig_genes_q20)))

  # Transcript-level DTU
  sig_tx_q05 <- res_tx[!is.na(res_tx$adj_pvalue) & res_tx$adj_pvalue < 0.05, ]
  sig_tx_q10 <- res_tx[!is.na(res_tx$adj_pvalue) & res_tx$adj_pvalue < 0.10, ]

  cat(sprintf("    Tx-level: q<0.05=%d, q<0.10=%d\n",
    nrow(sig_tx_q05), nrow(sig_tx_q10)))

  # Compute IF/dIF from TMM-normalized CPM
  dge_multi <- calcNormFactors(dge_multi, method = "TMM")
  cpm_vals <- cpm(dge_multi, log = FALSE)
  gene_totals_cpm <- rowsum(cpm_vals, group = gene_multi)
  if_values <- cpm_vals / gene_totals_cpm[gene_multi, ]
  if_values[is.nan(if_values)] <- 0

  idx_c1 <- which(samples_comp$timepoint == c1)
  idx_c2 <- which(samples_comp$timepoint == c2)

  mean_if_c1 <- rowMeans(if_values[, idx_c1, drop = FALSE])
  mean_if_c2 <- rowMeans(if_values[, idx_c2, drop = FALSE])
  dIF <- mean_if_c2 - mean_if_c1

  # Merge dIF with transcript results
  res_tx$mean_IF_c1 <- mean_if_c1[res_tx$feature_id]
  res_tx$mean_IF_c2 <- mean_if_c2[res_tx$feature_id]
  res_tx$dIF <- dIF[res_tx$feature_id]

  # Add gene name
  if (length(gene_to_name) > 0) {
    res_tx$gene_name <- gene_to_name[res_tx$gene_id]
    res_tx$gene_name[is.na(res_tx$gene_name)] <- res_tx$gene_id[is.na(res_tx$gene_name)]
    res_gene$gene_name <- gene_to_name[res_gene$gene_id]
    res_gene$gene_name[is.na(res_gene$gene_name)] <- res_gene$gene_id[is.na(res_gene$gene_name)]
  }

  # Combined DTU: gene q<0.05 AND |dIF|>=0.05
  sig_combined <- res_tx[res_tx$gene_id %in% sig_genes_q05$gene_id &
                         !is.na(res_tx$adj_pvalue) &
                         abs(res_tx$dIF) >= 0.05, ]
  sig_combined <- sig_combined[order(sig_combined$adj_pvalue), ]

  cat(sprintf("    Combined (gene q<0.05 + |dIF|>=0.05): %d tx in %d genes\n",
    nrow(sig_combined), length(unique(sig_combined$gene_id))))

  # Top DTU genes
  if (nrow(sig_genes_q05) > 0) {
    cat(sprintf("\n    Top DTU genes in %s:\n", cn))
    top_genes <- head(sig_genes_q05[order(sig_genes_q05$adj_pvalue), ], 15)
    for (j in 1:nrow(top_genes)) {
      gn <- ifelse(is.na(top_genes$gene_name[j]), top_genes$gene_id[j], top_genes$gene_name[j])
      cat(sprintf("    %2d. %-20s q=%.2e\n", j, gn, top_genes$adj_pvalue[j]))
    }
  }

  # Save
  comp_out <- file.path(OUTDIR, cn)
  dir.create(comp_out, showWarnings = FALSE)
  write.csv(res_gene, file.path(comp_out, "gene_results.csv"), row.names = FALSE)
  write.csv(res_tx, file.path(comp_out, "transcript_results.csv"), row.names = FALSE)
  if (nrow(sig_genes_q05) > 0)
    write.csv(sig_genes_q05, file.path(comp_out, "significant_genes_q005.csv"), row.names = FALSE)
  if (nrow(sig_combined) > 0)
    write.csv(sig_combined, file.path(comp_out, "significant_dtu_combined.csv"), row.names = FALSE)

  all_results[[cn]] <- list(
    comparison = cn,
    n_genes_tested = nrow(res_gene),
    n_tx_tested = nrow(res_tx),
    n_genes_q005 = nrow(sig_genes_q05),
    n_genes_q010 = nrow(sig_genes_q10),
    n_genes_q020 = nrow(sig_genes_q20),
    n_tx_q005 = nrow(sig_tx_q05),
    n_tx_q010 = nrow(sig_tx_q10),
    n_tx_combined = nrow(sig_combined)
  )

  rm(dge, d, dge_multi); gc(full = TRUE)
}

# ── 5. Summary ─────────────────────────────────────────────────────
cat("\n[5/5] Summary\n\n")
cat("  ──────────────────────────────────────────────────────────\n")
cat("  DRIMSeq DTU v6 SUMMARY (beta-binomial, per-comp filter)\n")
cat("  ──────────────────────────────────────────────────────────\n\n")
cat(sprintf("  %-25s %8s %8s %8s %8s\n",
  "Comparison", "Genes", "DTU_g05", "DTU_g10", "DTU_g20"))
cat("  "); cat(rep("-", 56), sep = ""); cat("\n")

overall <- data.frame()
for (cn in comp_df$name) {
  r <- all_results[[cn]]
  cat(sprintf("  %-25s %8d %8d %8d %8d\n",
    r$comparison, r$n_genes_tested,
    r$n_genes_q005, r$n_genes_q010, r$n_genes_q020))
  overall <- rbind(overall, data.frame(
    Comparison = r$comparison,
    Genes_tested = r$n_genes_tested,
    Tx_tested = r$n_tx_tested,
    Genes_q005 = r$n_genes_q005,
    Genes_q010 = r$n_genes_q010,
    Genes_q020 = r$n_genes_q020,
    Tx_q005 = r$n_tx_q005,
    Tx_q010 = r$n_tx_q010,
    Tx_combined = r$n_tx_combined,
    stringsAsFactors = FALSE
  ))
}
write.csv(overall, file.path(OUTDIR, "summary.csv"), row.names = FALSE)

# Cross-comparison
cat("\nCross-comparison:\n")
all_dtu_genes <- list()
for (cn in comp_df$name) {
  sig_file <- file.path(OUTDIR, cn, "significant_genes_q005.csv")
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
  cat(sprintf("  %d genes with DTU in >=1 comparison\n", nrow(gene_comp_count)))
  cat(sprintf("  %d genes with DTU in >=2 comparisons\n",
    sum(gene_comp_count$n_comparisons >= 2)))
  write.csv(gene_comp_count, file.path(OUTDIR, "multiswitch_genes.csv"), row.names = FALSE)
} else {
  cat("  No significant DTU at q<0.05\n")
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("  DRIMSeq v6 pipeline COMPLETE\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("════════════════════════════════════════════════════════════\n")

sink()