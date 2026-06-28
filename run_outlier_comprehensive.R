#!/usr/bin/env Rscript
# CTVT: Comprehensive outlier detection — multiple methods cross-validation
# After: Salmon + biotype + ~dog_id+timepoint + removeBatchEffect(dog_id)
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(tximport)
  library(pheatmap)
  library(matrixStats)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"

# Load data (same as rPCA)
tx2gene_pc <- read.csv(file.path(BASE, "tx2gene_protein_coding.csv"), stringsAsFactors = FALSE)
tx2gene_lnc <- read.csv(file.path(BASE, "tx2gene_lncrna_known.csv"), stringsAsFactors = FALSE)
tx2gene_biotype <- rbind(tx2gene_pc, tx2gene_lnc)

all_samples <- basename(list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE))
sample_meta <- data.frame(sampleID = all_samples, stringsAsFactors = FALSE)
sample_meta$timepoint <- factor(sapply(sample_meta$sampleID, function(s) {
  c(A="Day0",B="Day2",C="Day6",N="Recovered")[substr(s,3,3)]
}), levels = c("Day0", "Day2", "Day6", "Recovered"))
sample_meta$dog_id <- factor(sapply(sample_meta$sampleID, function(s) {
  substr(strsplit(s, "_")[[1]][1], 1, 2)
}))

files <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
names(files) <- sample_meta$sampleID
txi <- tximport(files, type = "salmon", tx2gene = tx2gene_biotype, countsFromAbundance = "lengthScaledTPM")
keep <- rowSums(txi$counts >= 10) >= 4
counts_raw <- round(txi$counts[keep, ])

dds <- DESeqDataSetFromMatrix(counts_raw, sample_meta, design = ~ dog_id + timepoint)
dds <- DESeq(dds, quiet = TRUE)
vsd <- vst(dds, blind = FALSE)
log_expr <- assay(vsd)
log_corr <- removeBatchEffect(log_expr, batch = sample_meta$dog_id)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  COMPREHENSIVE OUTLIER DETECTION — 5 METHODS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════════════
# METHOD 1: Sample-sample Euclidean distance (batch-corrected)
# ══════════════════════════════════════════════════════════════════
cat("── METHOD 1: Sample-sample distance (batch-corrected VST) ──\n")
sample_dist <- dist(t(log_corr))
dist_mat <- as.matrix(sample_dist)

# For each sample, mean distance to other samples of SAME timepoint
cat("\n  Mean distance to same-timepoint samples (lower = better):\n")
for (tp in c("Day0", "Day2", "Day6", "Recovered")) {
  tp_samples <- sample_meta$sampleID[sample_meta$timepoint == tp]
  tp_idx <- which(colnames(log_corr) %in% tp_samples)
  tp_dist <- dist_mat[tp_idx, tp_idx]
  diag(tp_dist) <- NA
  mean_dists <- colMeans(tp_dist, na.rm = TRUE)
  overall_mean <- mean(mean_dists)
  sd_mean <- sd(mean_dists)
  cat(sprintf("\n  %s (overall mean: %.1f, sd: %.1f):\n", tp, overall_mean, sd_mean))
  for (i in 1:length(mean_dists)) {
    s <- names(mean_dists)[i]
    z <- (mean_dists[i] - overall_mean) / sd_mean
    flag <- ifelse(abs(z) > 2, " ⚠️ HIGH", ifelse(abs(z) > 1.5, " ~moderate", ""))
    cat(sprintf("    %s: mean_dist=%.1f (z=%.2f)%s\n", s, mean_dists[i], z, flag))
  }
}

# ══════════════════════════════════════════════════════════════════
# METHOD 2: Correlation with timepoint group mean
# ══════════════════════════════════════════════════════════════════
cat("\n\n── METHOD 2: Correlation with timepoint group mean ──\n")
for (tp in c("Day0", "Day2", "Day6", "Recovered")) {
  tp_samples <- sample_meta$sampleID[sample_meta$timepoint == tp]
  tp_idx <- which(colnames(log_corr) %in% tp_samples)
  group_mean <- rowMeans(log_corr[, tp_idx])
  
  cat(sprintf("\n  %s:\n", tp))
  cors <- c()
  for (i in tp_idx) {
    s <- colnames(log_corr)[i]
    r <- cor(log_corr[, i], group_mean)
    cors <- c(cors, r)
    names(cors)[length(cors)] <- s
  }
  mean_r <- mean(cors)
  sd_r <- sd(cors)
  for (i in 1:length(cors)) {
    s <- names(cors)[i]
    z <- (cors[i] - mean_r) / sd_r
    flag <- ifelse(z < -2, " ⚠️ LOW CORR", ifelse(z < -1.5, " ~moderate", ""))
    cat(sprintf("    %s: r=%.4f (z=%.2f)%s\n", s, cors[i], z, flag))
  }
}

# ══════════════════════════════════════════════════════════════════
# METHOD 3: Hierarchical clustering (does sample cluster with its group?)
# ══════════════════════════════════════════════════════════════════
cat("\n\n── METHOD 3: Hierarchical clustering (batch-corrected) ──\n")
hc <- hclust(dist(t(log_corr)))
# For each sample, check if nearest neighbor is same timepoint
cat("  Nearest neighbor check (does sample cluster with same timepoint?):\n")
for (i in 1:ncol(log_corr)) {
  s <- colnames(log_corr)[i]
  tp_s <- sample_meta$timepoint[sample_meta$sampleID == s]
  # Find nearest neighbor
  dists <- dist_mat[s, ]
  dists[s] <- Inf
  nearest <- names(which.min(dists))
  tp_nearest <- sample_meta$timepoint[sample_meta$sampleID == nearest]
  same <- ifelse(tp_s == tp_nearest, "✓", "✗ WRONG GROUP")
  cat(sprintf("    %s (%s) → nearest: %s (%s) %s\n", s, tp_s, nearest, tp_nearest, same))
}

# ══════════════════════════════════════════════════════════════════
# METHOD 4: rPCA with 99% cutoff (stricter)
# ══════════════════════════════════════════════════════════════════
cat("\n\n── METHOD 4: rPCA with 99% cutoff (more conservative) ──\n")
library(rrcov)
rv <- rowVars(log_corr)
select_500 <- head(order(rv, decreasing = TRUE), 500)
pca_data <- t(log_corr[select_500, ])

rpca_99 <- PcaGrid(pca_data, k = 5, scale = FALSE, method = "mad", crit.pca.distances = 0.99)
od99 <- rpca_99$od
sd99 <- rpca_99$sd
cutoff_od99 <- rpca_99$cutoff.od
cutoff_sd99 <- rpca_99$cutoff.sd

cat(sprintf("  OD cutoff (99%%): %.2f | SD cutoff (99%%): %.2f\n", cutoff_od99, cutoff_sd99))

results99 <- data.frame(
  Sample = colnames(log_corr),
  Timepoint = sample_meta$timepoint,
  Dog = sample_meta$dog_id,
  OD99 = round(od99, 2),
  OD99_ratio = round(od99 / cutoff_od99, 2),
  SD99 = round(sd99, 2),
  SD99_ratio = round(sd99 / cutoff_sd99, 2),
  stringsAsFactors = FALSE
)
results99 <- results99[order(-pmax(results99$OD99_ratio, results99$SD99_ratio)), ]

flagged_99 <- results99[results99$OD99_ratio > 1 | results99$SD99_ratio > 1, ]
cat(sprintf("  Flagged at 99%%: %d samples\n", nrow(flagged_99)))
if (nrow(flagged_99) > 0) {
  for (i in 1:min(5, nrow(flagged_99))) {
    r <- flagged_99[i, ]
    cat(sprintf("    %s Dog=%s Time=%s  OD=%.2fx  SD=%.2fx\n", r$Sample, r$Dog, r$Timepoint, r$OD99_ratio, r$SD99_ratio))
  }
}

# ══════════════════════════════════════════════════════════════════
# METHOD 5: Cook's distance per sample (DESeq2 built-in)
# ══════════════════════════════════════════════════════════════════
cat("\n\n── METHOD 5: Cook's distance per sample (DESeq2) ──\n")
cooks <- cooksDistance(dds)
# Sum of Cook's distance per sample (high = sample has many outlier genes)
sample_cooks <- colSums(cooks > 0.5, na.rm = TRUE)
sample_cooks_total <- colSums(cooks, na.rm = TRUE)
results_cook <- data.frame(
  Sample = names(sample_cooks),
  Timepoint = sample_meta$timepoint,
  Dog = sample_meta$dog_id,
  Cooks_total = round(sample_cooks_total, 1),
  Cooks_flagged_genes = sample_cooks,
  stringsAsFactors = FALSE
)
results_cook <- results_cook[order(-results_cook$Cooks_flagged_genes), ]
cat("  Genes with Cook's distance > 0.5 per sample (top 10):\n")
for (i in 1:min(10, nrow(results_cook))) {
  r <- results_cook[i, ]
  cat(sprintf("    %s Dog=%s Time=%s  flagged_genes=%d  total_cooks=%.1f\n", 
      r$Sample, r$Dog, r$Timepoint, r$Cooks_flagged_genes, r$Cooks_total))
}

# ══════════════════════════════════════════════════════════════════
# SUMMARY: Consensus across all methods
# ══════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  CONSENSUS SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Flag samples based on consensus
# A sample is "concerning" if flagged by 3+ methods
# A sample is "borderline" if flagged by 2 methods
# A sample is "fine" if flagged by 0-1 methods

# Save distance heatmap
pdf(file.path(BASE, "sample_distance_heatmap.pdf"), width = 12, height = 10)
annotation_row <- data.frame(
  Timepoint = sample_meta$timepoint,
  Dog = sample_meta$dog_id,
  row.names = sample_meta$sampleID
)
pheatmap(dist_mat, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 7, fontsize_col = 7,
         main = "Sample-sample distance (batch-corrected VST)")
dev.off()

cat("Saved: sample_distance_heatmap.pdf\n")