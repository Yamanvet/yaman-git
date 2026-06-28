#!/usr/bin/env Rscript
# ======================================================================
# CTVT: DTU Pipeline v7 — limma diffSplice + duplicateCorrelation
#
# Based on deep methodological review (see analysis PDF):
#   - diffSplice tests DTU correctly (transcript logFC vs gene average)
#   - duplicateCorrelation models dog_id as RANDOM effect (preserves df)
#   - Empirical Bayes variance moderation stabilizes small-n estimates
#   - Standard B-H FDR (no empirical null overcorrection)
#   - IF/dIF from TPM (length-adjusted), NOT from CPM (length-biased)
#   - Per-comparison filtering with relaxed thresholds (min.count=5, n.small=4)
#
# Note: catchSalmon requires Gibbs samples (--numGibbsSamples) which we
# don't have. Using tximport lengthScaledTPM for voom instead.
# The diffSplice + duplicateCorrelation approach still correctly tests
# DTU even without divided counts.
# ======================================================================

.libPaths(c("~/R/library", .libPaths()))

Sys.setenv(R_MAX_VSIZE = "500e9")
options(future.globals.maxSize = 10 * 1024^3)
gc(full = TRUE)

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(tximport)
  library(dplyr)
})

# ── Config ──────────────────────────────────────────────────────────
BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"
GTF_FILE <- file.path(BASE, "stringtie_v4/merged_stringtie_v4.gtf")
OUTDIR <- file.path(BASE, "diffsplice_v7")
OUTLIERS <- c("T1B_S2", "T8C_S10", "T9N_S26")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUTDIR, "run.log"), split = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CTVT: DTU v7 — limma diffSplice + duplicateCorrelation\n")
cat("  R:", as.character(getRversion()), "\n")
cat("  limma:", as.character(packageVersion("limma")), "\n")
cat("  edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── 1. Sample metadata ─────────────────────────────────────────────
cat("[1/7] Building sample metadata...\n")

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

# ── 2. Import via tximport ─────────────────────────────────────────
cat("\n[2/7] Importing Salmon quant via tximport...\n")

# tximport with lengthScaledTPM for voom (correct for DTU via diffSplice)
# and TPM for IF calculation (length-adjusted, as per analysis PDF)
files <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
names(files) <- sample_meta$sampleID

txi <- tximport(files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "lengthScaledTPM")

cat(sprintf("  Imported: %d transcripts x %d samples\n",
  nrow(txi$counts), ncol(txi$counts)))
cat(sprintf("  Abundance matrix: %d x %d\n", nrow(txi$abundance), ncol(txi$abundance)))

# Save raw imports
saveRDS(txi, file.path(OUTDIR, "01_tximport.rds"))
gc(full = TRUE)

# ── 3. GTF gene-transcript mapping ─────────────────────────────────
cat("\n[3/7] Importing GTF for gene-transcript mapping...\n")

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

# Build tx2gene for tximport gene-level summarization
tx2gene <- gene_tx_map[, c("transcript_id", "gene_id")]
tx2gene <- tx2gene[!duplicated(tx2gene), ]

# ── 4. Per-comparison DTU via diffSplice + duplicateCorrelation ─────
cat("\n[4/7] Per-comparison DTU (diffSplice + duplicateCorrelation)...\n")

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

  # Subset samples for this comparison
  keep_comp <- sample_meta$timepoint %in% c(c1, c2)
  samples_comp <- sample_meta[keep_comp, ]

  n1 <- sum(samples_comp$timepoint == c1)
  n2 <- sum(samples_comp$timepoint == c2)

  cat(sprintf("    %s (n=%d) vs %s (n=%d)\n", c1, n1, c2, n2))

  # Subset tximport data
  # txi$counts = lengthScaledTPM (since countsFromAbundance="lengthScaledTPM")
  # txi$abundance = TPM (for IF/dIF calculation)
  keep_cols <- which(colnames(txi$counts) %in% samples_comp$sampleID)
  lstp_comp <- txi$counts[, keep_cols, drop = FALSE]  # lengthScaledTPM for voom
  abundance_comp <- txi$abundance[, keep_cols, drop = FALSE]  # TPM for IF

  # Per-comparison filter: relaxed thresholds (min.count=5, n.small=4)
  # Only keep transcripts expressed in at least min_n samples
  n_small <- min(n1, n2)
  min_n <- max(4, floor(n_small / 2))  # at least 4 or half the smaller group

  # Filter on lengthScaledTPM counts (voom input)
  keep_idx <- rowSums(lstp_comp >= 5) >= min_n
  lstp_filt <- lstp_comp[keep_idx, ]
  abundance_filt <- abundance_comp[keep_idx, ]

  cat(sprintf("    After relaxed filter (lstp>=5 in %d/%d samples): %d tx\n",
    min_n, ncol(lstp_filt), nrow(lstp_filt)))

  # Map to genes
  tx_names <- rownames(lstp_filt)
  gene_ids <- tx_to_gene[tx_names]
  unmapped <- is.na(gene_ids)
  gene_ids[unmapped] <- tx_names[unmapped]

  # Keep only multi-isoform genes
  gene_counts <- table(gene_ids)
  multi_genes <- names(gene_counts[gene_counts >= 2])
  keep_multi <- gene_ids %in% multi_genes

  lstp_multi <- lstp_filt[keep_multi, ]
  abundance_multi <- abundance_filt[keep_multi, ]
  gene_multi <- gene_ids[keep_multi]
  tx_multi <- tx_names[keep_multi]

  cat(sprintf("    Multi-isoform: %d genes (%d tx)\n",
    length(multi_genes), sum(keep_multi)))

  if (length(multi_genes) < 10) {
    cat("    Too few multi-isoform genes, skipping\n")
    next
  }

  # ── DGEList from lengthScaledTPM (for voom) ──
  # lstp = lengthScaledTPM counts (length-adjusted, suitable for voom)
  dge <- DGEList(counts = lstp_multi)
  dge <- calcNormFactors(dge, method = "TMM")

  # ── Design matrix (condition only) ──
  # duplicateCorrelation will handle dog_id as random effect
  samples_comp$condition <- factor(samples_comp$timepoint, levels = c(c1, c2))
  design <- model.matrix(~ condition, data = samples_comp)
  colnames(design) <- c("Intercept", paste0("condition", c2))

  # ── voom with quality weights ──
  cat("    voomWithQualityWeights...\n")
  v <- voomWithQualityWeights(dge, design, plot = FALSE)

  # ── duplicateCorrelation for dog_id (random effect) ──
  cat("    duplicateCorrelation (dog_id as random effect)...\n")
  corfit <- duplicateCorrelation(v$E, design, block = samples_comp$dog_id)
  cat(sprintf("    Intra-dog correlation: %.4f\n", corfit$consensus.correlation))

  # ── Re-fit voom with correlation ──
  cat("    voom with consensus correlation...\n")
  v <- voomWithQualityWeights(dge, design, block = samples_comp$dog_id,
                               correlation = corfit$consensus.correlation, plot = FALSE)

  # ── Fit linear model ──
  cat("    lmFit...\n")
  fit <- lmFit(v, design, block = samples_comp$dog_id,
               correlation = corfit$consensus.correlation)

  # ── diffSplice ──
  cat("    diffSplice...\n")
  # diffSplice: verbose=TRUE shows progress; geneid/transcriptid match tx to genes
  sp <- diffSplice(fit, geneid = gene_multi, exonid = tx_multi,
                    robust = FALSE, verbose = TRUE)

  # ── Extract results ──
  cat("    topSplice (exon-level = transcript-level)...\n")
  tx_results <- topSplice(sp, coef = paste0("condition", c2), test = "t",
                            number = Inf, sort.by = "none")

  # Gene-level Simes test (lowercase!)
  gene_results <- topSplice(sp, coef = paste0("condition", c2), test = "simes",
                             number = Inf, sort.by = "none")

  cat(sprintf("    Transcripts tested: %d\n", nrow(tx_results)))
  cat(sprintf("    Genes tested: %d\n", length(unique(gene_results$GeneID))))

  # ── IF/dIF from TPM (length-adjusted, NOT CPM) ──
  # Per the analysis PDF: IF must be calculated from TPM/length-adjusted values
  # not from TMM-CPM which has length bias
  gene_totals_tpm <- rowsum(abundance_multi, group = gene_multi)
  if_values <- abundance_multi / gene_totals_tpm[gene_multi, ]
  if_values[is.nan(if_values)] <- 0

  idx_c1 <- which(samples_comp$timepoint == c1)
  idx_c2 <- which(samples_comp$timepoint == c2)

  mean_if_c1 <- rowMeans(if_values[, idx_c1, drop = FALSE])
  mean_if_c2 <- rowMeans(if_values[, idx_c2, drop = FALSE])
  dIF <- mean_if_c2 - mean_if_c1

  # ── Merge dIF into results ──
  # diffSplice uses exon/transcript IDs as rownames
  # dIF indexed by transcript ID (ExonID column matches our tx_multi)
  names(dIF) <- tx_multi
  names(mean_if_c1) <- tx_multi
  names(mean_if_c2) <- tx_multi

  tx_results$dIF <- dIF[tx_results$ExonID]
  tx_results$mean_IF_c1 <- mean_if_c1[tx_results$ExonID]
  tx_results$mean_IF_c2 <- mean_if_c2[tx_results$ExonID]
  tx_results$gene_name <- gene_to_name[tx_results$GeneID]
  tx_results$gene_name[is.na(tx_results$gene_name)] <- tx_results$GeneID[is.na(tx_results$gene_name)]

  # Gene results
  gene_results$gene_name <- gene_to_name[gene_results$GeneID]
  gene_results$gene_name[is.na(gene_results$gene_name)] <- gene_results$GeneID[is.na(gene_results$gene_name)]

  # ── Significance at multiple thresholds ──
  for (thresh in c(0.05, 0.10, 0.20)) {
    # Transcript-level (t-test)
    sig_tx <- tx_results[!is.na(tx_results$FDR) &
                          tx_results$FDR < thresh &
                          abs(tx_results$dIF) >= 0.05, ]
    # Gene-level (Simes)
    sig_genes <- gene_results[!is.na(gene_results$FDR) &
                               gene_results$FDR < thresh, ]

    cat(sprintf("    q<%.2f + |dIF|>=0.05: %d tx / %d genes (Simes: %d genes)\n",
      thresh, nrow(sig_tx), length(unique(sig_tx$GeneID)), nrow(sig_genes)))
  }

  # Primary results at q<0.05 + |dIF|>=0.05
  sig_q05 <- tx_results[!is.na(tx_results$FDR) &
                         tx_results$FDR < 0.05 &
                         abs(tx_results$dIF) >= 0.05, ]
  sig_q05 <- sig_q05[order(sig_q05$FDR), ]

  # Top DTU transcripts
  if (nrow(sig_q05) > 0) {
    cat(sprintf("\n    Top DTU in %s:\n", cn))
    top <- head(sig_q05, 15)
    for (j in 1:nrow(top)) {
      gn <- ifelse(is.na(top$gene_name[j]), top$GeneID[j], top$gene_name[j])
      cat(sprintf("    %2d. %-20s iso=%-20s dIF=%+.4f q=%.2e\n",
        j, gn, top$ExonID[j], top$dIF[j], top$FDR[j]))
    }
  }

  # Save
  comp_out <- file.path(OUTDIR, cn)
  dir.create(comp_out, showWarnings = FALSE)
  write.csv(tx_results, file.path(comp_out, "transcript_results.csv"), row.names = TRUE)
  write.csv(gene_results, file.path(comp_out, "gene_results.csv"), row.names = TRUE)
  if (nrow(sig_q05) > 0)
    write.csv(sig_q05, file.path(comp_out, "significant_dtu_q005.csv"), row.names = TRUE)

  # Simes gene-level
  sig_genes_q05 <- gene_results[!is.na(gene_results$FDR) & gene_results$FDR < 0.05, ]
  if (nrow(sig_genes_q05) > 0)
    write.csv(sig_genes_q05, file.path(comp_out, "significant_genes_q005.csv"), row.names = TRUE)

  all_results[[cn]] <- list(
    comparison = cn,
    n_tx_tested = nrow(tx_results),
    n_genes_tested = length(unique(gene_results$Gene)),
    consensus_cor = corfit$consensus.correlation,
    n_tx_q005 = sum(!is.na(tx_results$FDR) & tx_results$FDR < 0.05 & abs(tx_results$dIF) >= 0.05),
    n_genes_q005_simes = sum(!is.na(gene_results$FDR) & gene_results$FDR < 0.05),
    n_tx_q010 = sum(!is.na(tx_results$FDR) & tx_results$FDR < 0.10 & abs(tx_results$dIF) >= 0.05),
    n_tx_q020 = sum(!is.na(tx_results$FDR) & tx_results$FDR < 0.20 & abs(tx_results$dIF) >= 0.05)
  )

  rm(dge, v, corfit, fit, sp); gc(full = TRUE)
}

# ── 5. Summary ─────────────────────────────────────────────────────
cat("\n[5/7] Summary\n\n")
cat("  ──────────────────────────────────────────────────────────\n")
cat("  DTU v7 SUMMARY (diffSplice + duplicateCorrelation)\n")
cat("  ──────────────────────────────────────────────────────────\n\n")
cat(sprintf("  %-25s %8s %8s %8s %8s\n",
  "Comparison", "Tx_test", "DTU(q05)", "Genes(Simes)", "DTU(q10)"))
cat("  "); cat(rep("-", 62), sep = ""); cat("\n")

overall <- data.frame()
for (cn in comp_df$name) {
  r <- all_results[[cn]]
  cat(sprintf("  %-25s %8d %8d %8d %8d\n",
    r$comparison, r$n_tx_tested, r$n_tx_q005,
    r$n_genes_q005_simes, r$n_tx_q010))
  overall <- rbind(overall, data.frame(
    Comparison = r$comparison,
    Tx_tested = r$n_tx_tested,
    Genes_tested = r$n_genes_tested,
    Consensus_cor = r$consensus_cor,
    DTU_q005 = r$n_tx_q005,
    Genes_Simes_q005 = r$n_genes_q005_simes,
    DTU_q010 = r$n_tx_q010,
    DTU_q020 = r$n_tx_q020,
    stringsAsFactors = FALSE
  ))
}
write.csv(overall, file.path(OUTDIR, "summary.csv"), row.names = FALSE)

# ── 6. Cross-comparison ────────────────────────────────────────────
cat("\n[6/7] Cross-comparison: genes with DTU in multiple contrasts\n")

all_dtu_genes <- list()
for (cn in comp_df$name) {
  sig_file <- file.path(OUTDIR, cn, "significant_genes_q005.csv")
  if (file.exists(sig_file)) {
    df <- read.csv(sig_file)
    if (nrow(df) > 0) all_dtu_genes[[cn]] <- unique(df$Gene)
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
  if (sum(gene_comp_count$n_comparisons >= 3) > 0)
    cat(sprintf("  %d genes with DTU in >=3 comparisons\n",
      sum(gene_comp_count$n_comparisons >= 3)))
  write.csv(gene_comp_count, file.path(OUTDIR, "multiswitch_genes.csv"), row.names = FALSE)
} else {
  cat("  No significant DTU at q<0.05\n")
}

# ── 7. Gene-level pathway annotation ────────────────────────────────
cat("\n[7/7] Gene annotation summary\n")

all_sig <- list()
for (cn in comp_df$name) {
  sig_file <- file.path(OUTDIR, cn, "significant_dtu_q005.csv")
  if (file.exists(sig_file)) {
    df <- read.csv(sig_file)
    if (nrow(df) > 0) {
      df$comparison <- cn
      all_sig[[cn]] <- df
    }
  }
}

if (length(all_sig) > 0) {
  combined <- do.call(rbind, all_sig)
  cat(sprintf("  Total significant DTU events: %d\n", nrow(combined)))
  cat(sprintf("  Unique genes: %d\n", length(unique(combined$GeneID))))
  cat(sprintf("  Unique transcripts: %d\n", length(unique(combined$ExonID))))

  # Per-comparison breakdown
  for (cn in comp_df$name) {
    if (cn %in% names(all_sig)) {
      df <- all_sig[[cn]]
      cat(sprintf("\n  %s: %d DTU in %d genes\n", cn, nrow(df), length(unique(df$GeneID))))
      if (nrow(df) <= 20) {
        for (j in 1:nrow(df)) {
          gn <- ifelse(is.na(df$gene_name[j]), df$GeneID[j], df$gene_name[j])
          cat(sprintf("    %2d. %-20s dIF=%+.4f q=%.2e\n",
            j, gn, df$dIF[j], df$FDR[j]))
        }
      }
    }
  }
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("  DTU v7 pipeline COMPLETE\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("  Key: diffSplice tests DTU (transcript vs gene logFC)\n")
cat("  Key: duplicateCorrelation models dog_id as random effect\n")
cat("  Key: IF/dIF from TPM (length-adjusted, NOT CPM)\n")
cat("════════════════════════════════════════════════════════════\n")

sink()