#!/usr/bin/env Rscript
# ======================================================================
# CTVT: satuRn DTU v4 — Correct satuRn API usage
#
# KEY FIX: Previous v2/v3 used voom+limma on counts → tests DGE, NOT DTU!
# satuRn uses quasi-binomial GLM on isoform PROPORTIONS within each gene.
#
# satuRn results structure:
#   fitDTU → adds rowData$fitDTUModels (StatModel objects)
#   testDTU → adds rowData$fitDTUResult_<contrastName> with columns:
#     estimates, se, df, t, pval, regular_FDR, empirical_pval, empirical_FDR
#   estimates = log-odds ratio of isoform usage
# ======================================================================

.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_MAX_VSIZE = "500e9")
options(future.globals.maxSize = 10 * 1024^3)
gc(full = TRUE)

suppressPackageStartupMessages({
  library(satuRn)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(edgeR)
})

# ── Config ──────────────────────────────────────────────────────────
BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"
GTF_FILE <- file.path(BASE, "stringtie_v4/merged_stringtie_v4.gtf")
OUTDIR <- file.path(BASE, "saturn_v4")
OUTLIERS <- c("T1B_S2", "T8C_S10", "T9N_S26")

ncores <- min(parallel::detectCores(), 60)
BPPARAM <- MulticoreParam(workers = ncores, progressbar = TRUE)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUTDIR, "run.log"), split = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CTVT: satuRn DTU v4 — quasi-binomial GLM on proportions\n")
cat("  R:", as.character(getRversion()), "\n")
cat("  satuRn:", as.character(packageVersion("satuRn")), "\n")
cat("  GTF:", GTF_FILE, "\n")
cat("  Quant:", QUANT_DIR, "\n")
cat("  Method: fitDTU (quasi-binomial GLM) → testDTU (contrast)\n")
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

# Raw NumReads matrix
raw_counts <- matrix(0, nrow = n_tx, ncol = nrow(sample_meta),
                      dimnames = list(tx_names, sample_meta$sampleID))

for (i in 1:nrow(sample_meta)) {
  df <- read.delim(sample_meta$path[i], stringsAsFactors = FALSE)
  raw_counts[df$Name, i] <- df$NumReads
  if (i %% 10 == 0) cat(sprintf("    Read %d/%d samples\n", i, nrow(sample_meta)))
}
cat(sprintf("    Read %d/%d — done\n", i, nrow(sample_meta)))

# Gene mapping for all transcripts
gene_for_tx <- tx_to_gene[tx_names]
unmapped <- is.na(gene_for_tx)
gene_for_tx[unmapped] <- tx_names[unmapped]

saveRDS(list(raw_counts = raw_counts, gene_for_tx = gene_for_tx),
        file.path(OUTDIR, "01_salmon_counts.rds"))
cat("  Saved: 01_salmon_counts.rds\n")
gc(full = TRUE)

# ── 4. Per-comparison satuRn DTU ───────────────────────────────────
cat("\n[4/6] Per-comparison satuRn DTU (quasi-binomial GLM)...\n")

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

  # Subset samples
  keep_comp <- sample_meta$timepoint %in% c(c1, c2)
  samples_comp <- sample_meta[keep_comp, ]

  cat(sprintf("    %s (n=%d) vs %s (n=%d)\n",
    c1, sum(samples_comp$timepoint == c1),
    c2, sum(samples_comp$timepoint == c2)))

  # Subset raw counts
  raw_comp <- raw_counts[, samples_comp$sampleID]

  # DGEList for filtering
  dge <- DGEList(counts = raw_comp)
  dge$samples$group <- factor(samples_comp$timepoint, levels = c(c1, c2))

  # Filter per-comparison
  keep <- filterByExpr(dge, group = dge$samples$group, min.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  filtered_tx <- rownames(dge$counts)
  gene_ids_filt <- gene_for_tx[filtered_tx]

  n_retained <- length(filtered_tx)
  cat(sprintf("    After per-comparison filter: %d / %d tx\n", n_retained, n_tx))

  # Keep only multi-isoform genes (DTU requires >= 2 isoforms per gene)
  gene_counts <- table(gene_ids_filt)
  multi_genes <- names(gene_counts[gene_counts >= 2])
  keep_multi <- gene_ids_filt %in% multi_genes

  n_multi <- sum(keep_multi)
  cat(sprintf("    Multi-isoform genes: %d genes (%d tx)\n",
    length(multi_genes), n_multi))

  if (n_multi == 0) {
    cat("    No multi-isoform genes — skipping\n")
    all_results[[cn]] <- list(comparison = cn, n_retained = n_retained,
      n_multi_genes = 0, n_multi_tx = 0, n_dtu_q005 = 0, n_genes_q005 = 0,
      n_dtu_q010 = 0, n_genes_q010 = 0)
    next
  }

  # Subset to multi-isoform
  dge_multi <- dge[keep_multi, ]
  tx_multi <- rownames(dge_multi$counts)
  gene_multi <- gene_for_tx[tx_multi]

  # ── Create SummarizedExperiment for satuRn ──────────────────────
  col_data <- DataFrame(
    condition = factor(samples_comp$timepoint, levels = c(c1, c2)),
    dog_id = samples_comp$dog_id,
    row.names = samples_comp$sampleID
  )

  row_data <- DataFrame(
    isoform_id = tx_multi,
    gene_id = gene_multi,
    row.names = tx_multi
  )

  se <- SummarizedExperiment(
    assays = list(counts = dge_multi$counts),
    rowData = row_data,
    colData = col_data
  )

  cat(sprintf("    SE: %d tx x %d samples\n", nrow(se), ncol(se)))

  # ── fitDTU ────────────────────────────────────────────────────────
  cat("    fitDTU (quasi-binomial GLM)...\n")
  fit_time <- system.time({
    se_fit <- fitDTU(se, formula = ~ dog_id + condition,
                     parallel = TRUE, BPPARAM = BPPARAM)
  })
  cat(sprintf("    fitDTU took %.1f sec\n", fit_time["elapsed"]))

  # Save fit object
  saveRDS(se_fit, file.path(OUTDIR, paste0(cn, "_fit.rds")))

  # ── Build contrast matrix ─────────────────────────────────────────
  # Get coefficient names from the fitted model
  first_model <- rowData(se_fit)$fitDTUModels[[1]]
  coef_names <- names(first_model@params$coefficients)
  cat(sprintf("    Coefficients: %s\n", paste(coef_names, collapse = ", ")))

  # satuRn prefixes coefficients with design e.g. designconditionDay2
  contrast_coef <- paste0("designcondition", c2)

  if (!(contrast_coef %in% coef_names)) {
    cat(sprintf("    WARNING: %s not in coefficients!\n", contrast_coef))
    cat("    Available coefficients:", paste(coef_names, collapse = ", "), "\n")
    all_results[[cn]] <- list(comparison = cn, n_retained = n_retained,
      n_multi_genes = length(multi_genes), n_multi_tx = n_multi,
      n_dtu_q005 = NA, n_genes_q005 = NA, n_dtu_q010 = NA, n_genes_q010 = NA,
      error = "Contrast coefficient not found")
    next
  }

  L <- matrix(0, nrow = length(coef_names), ncol = 1,
              dimnames = list(coef_names, cn))
  L[contrast_coef, 1] <- 1
  cat(sprintf("    Contrast: %s = 1 (baseline = %s)\n", contrast_coef, c1))

  # ── testDTU ───────────────────────────────────────────────────────
  cat("    testDTU...\n")
  test_time <- system.time({
    se_test <- testDTU(se_fit, contrasts = L, diagplot1 = FALSE,
                       diagplot2 = FALSE, sort = FALSE)
  })
  cat(sprintf("    testDTU took %.1f sec\n", test_time["elapsed"]))

  saveRDS(se_test, file.path(OUTDIR, paste0(cn, "_test.rds")))

  # ── Extract results ──────────────────────────────────────────────
  # Results are in rowData column named "fitDTUResult_<contrast_name>"
  result_col <- paste0("fitDTUResult_", cn)

  if (!(result_col %in% colnames(rowData(se_test)))) {
    cat(sprintf("    ERROR: Result column '%s' not found\n", result_col))
    cat("    Available columns:", paste(colnames(rowData(se_test)), collapse = ", "), "\n")
    all_results[[cn]] <- list(comparison = cn, n_retained = n_retained,
      n_multi_genes = length(multi_genes), n_multi_tx = n_multi,
      n_dtu_q005 = NA, n_genes_q005 = NA, error = "Result column not found")
    next
  }

  res_df <- as.data.frame(rowData(se_test)[[result_col]])
  res_df$transcript_id <- rownames(res_df)
  res_df$gene_id <- gene_multi

  if (length(gene_to_name) > 0) {
    res_df$gene_name <- gene_to_name[res_df$gene_id]
    res_df$gene_name[is.na(res_df$gene_name)] <- res_df$gene_id[is.na(res_df$gene_name)]
  }

  # ── Compute IF/dIF from TMM-normalized CPM ───────────────────────
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

  res_df$mean_IF_c1 <- mean_if_c1[res_df$transcript_id]
  res_df$mean_IF_c2 <- mean_if_c2[res_df$transcript_id]
  res_df$dIF <- dIF[res_df$transcript_id]

  # ── Significance ──────────────────────────────────────────────────
  # Use empirical_FDR if available (>500 features), else regular_FDR
  fdr_col <- if ("empirical_FDR" %in% colnames(res_df) &&
                   sum(!is.na(res_df$empirical_FDR)) > 100) {
    cat("    Using empirical FDR (empirical correction applied)\n")
    "empirical_FDR"
  } else {
    cat("    Using regular FDR (no empirical correction)\n")
    "regular_FDR"
  }

  sig_q05 <- res_df[!is.na(res_df[[fdr_col]]) &
                     res_df[[fdr_col]] < 0.05 &
                     abs(res_df$dIF) >= 0.05, ]
  sig_q05 <- sig_q05[order(sig_q05[[fdr_col]]), ]

  sig_q10 <- res_df[!is.na(res_df[[fdr_col]]) &
                     res_df[[fdr_col]] < 0.10 &
                     abs(res_df$dIF) >= 0.05, ]

  cat(sprintf("    Tested: %d tx\n", sum(!is.na(res_df$pval))))
  cat(sprintf("    Significant (q<0.05, |dIF|>=0.05): %d tx / %d genes\n",
    nrow(sig_q05), length(unique(sig_q05$gene_id))))
  cat(sprintf("    Relaxed   (q<0.10, |dIF|>=0.05): %d tx / %d genes\n",
    nrow(sig_q10), length(unique(sig_q10$gene_id))))

  # Also report without |dIF| threshold
  sig_q05_nodif <- res_df[!is.na(res_df[[fdr_col]]) &
                           res_df[[fdr_col]] < 0.05, ]
  cat(sprintf("    (without dIF filter): %d tx / %d genes\n",
    nrow(sig_q05_nodif), length(unique(sig_q05_nodif$gene_id))))

  if (nrow(sig_q05) > 0) {
    cat(sprintf("\n    Top DTU in %s:\n", cn))
    top <- head(sig_q05, 10)
    for (j in 1:nrow(top)) {
      gn <- ifelse(is.na(top$gene_name[j]), top$gene_id[j], top$gene_name[j])
      cat(sprintf("    %2d. %-18s iso=%-18s est=%+.4f dIF=%+.4f q=%.2e\n",
        j, gn, top$transcript_id[j], top$estimates[j], top$dIF[j], top[[fdr_col]][j]))
    }
  }

  # Save
  comp_out <- file.path(OUTDIR, cn)
  dir.create(comp_out, showWarnings = FALSE)
  write.csv(res_df, file.path(comp_out, "all_transcript_results.csv"), row.names = FALSE)
  if (nrow(sig_q05) > 0)
    write.csv(sig_q05, file.path(comp_out, "significant_dtu_q005.csv"), row.names = FALSE)
  if (nrow(sig_q10) > 0)
    write.csv(sig_q10, file.path(comp_out, "significant_dtu_q010.csv"), row.names = FALSE)

  all_results[[cn]] <- list(
    comparison = cn, n_retained = n_retained,
    n_multi_genes = length(multi_genes), n_multi_tx = n_multi,
    n_dtu_q005 = nrow(sig_q05), n_genes_q005 = length(unique(sig_q05$gene_id)),
    n_dtu_q010 = nrow(sig_q10), n_genes_q010 = length(unique(sig_q10$gene_id)),
    n_dtu_q005_nodif = nrow(sig_q05_nodif),
    fdr_used = fdr_col
  )

  rm(dge, se, se_fit, se_test, dge_multi); gc(full = TRUE)
}

# ── 5. Summary ─────────────────────────────────────────────────────
cat("\n[5/6] Summary\n\n")
cat("  ──────────────────────────────────────────────────────────\n")
cat("  satuRn DTU v4 SUMMARY (quasi-binomial GLM, per-comp filter)\n")
cat("  ──────────────────────────────────────────────────────────\n\n")
cat(sprintf("  %-25s %8s %8s %8s %8s\n",
  "Comparison", "MultiTx", "DTU(q05)", "Genes(q05)", "DTU(q10)"))
cat("  "); cat(rep("─", 56), sep = ""); cat("\n")

overall <- data.frame()
for (cn in comp_df$name) {
  r <- all_results[[cn]]
  cat(sprintf("  %-25s %8d %8d %8d %8d\n",
    r$comparison, r$n_multi_tx,
    r$n_dtu_q005, r$n_genes_q005, r$n_dtu_q010))
  overall <- rbind(overall, data.frame(
    Comparison = r$comparison, Retained = r$n_retained,
    MultiTx = r$n_multi_tx, MultiGenes = r$n_multi_genes,
    DTU_q005 = r$n_dtu_q005, Genes_q005 = r$n_genes_q005,
    DTU_q010 = r$n_dtu_q010, Genes_q010 = r$n_genes_q010,
    FDR_used = ifelse(is.null(r$fdr_used), NA, r$fdr_used),
    stringsAsFactors = FALSE
  ))
}
write.csv(overall, file.path(OUTDIR, "summary.csv"), row.names = FALSE)

# ── 6. Cross-comparison ────────────────────────────────────────────
cat("\n[6/6] Cross-comparison\n")

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
  cat(sprintf("  %d genes with DTU in >=1 comparison\n", nrow(gene_comp_count)))
  cat(sprintf("  %d genes with DTU in >=2 comparisons\n",
    sum(gene_comp_count$n_comparisons >= 2)))
  write.csv(gene_comp_count, file.path(OUTDIR, "multiswitch_genes.csv"),
    row.names = FALSE)
} else {
  cat("  No significant DTU at q<0.05\n")
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("  satuRn v4 pipeline COMPLETE\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("════════════════════════════════════════════════════════════\n")

sink()