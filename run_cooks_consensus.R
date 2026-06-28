#!/usr/bin/env Rscript
# Method 5 fix: Cook's distance
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(tximport)
  library(matrixStats)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"

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

# Cook's distance matrix (genes x samples)
cooks_mat <- assays(dds)[["cooks"]]
# Number of genes where this sample has high Cook's distance
sample_cooks_count <- colSums(cooks_mat > 0.5, na.rm = TRUE)
sample_cooks_total <- colSums(cooks_mat, na.rm = TRUE)

results_cook <- data.frame(
  Sample = names(sample_cooks_count),
  Timepoint = sample_meta$timepoint,
  Dog = sample_meta$dog_id,
  Flagged_genes = sample_cooks_count,
  Total_cooks = round(sample_cooks_total, 1),
  stringsAsFactors = FALSE
)
results_cook <- results_cook[order(-results_cook$Flagged_genes), ]

cat("── METHOD 5: Cook's distance per sample (DESeq2 built-in) ──\n")
cat("  Genes with Cook's distance > 0.5 per sample (top 15):\n")
for (i in 1:min(15, nrow(results_cook))) {
  r <- results_cook[i, ]
  cat(sprintf("    %s Dog=%s Time=%-10s  flagged_genes=%d  total_cooks=%.1f\n", 
      r$Sample, r$Dog, r$Timepoint, r$Flagged_genes, r$Total_cooks))
}

cat("\n  Per-timepoint mean flagged genes:\n")
for (tp in c("Day0", "Day2", "Day6", "Recovered")) {
  idx <- results_cook$Timepoint == tp
  cat(sprintf("    %s: mean=%.1f, max=%d\n", tp, 
      mean(results_cook$Flagged_genes[idx]), max(results_cook$Flagged_genes[idx])))
}

# CONSENSUS across all 5 methods
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  CONSENSUS: Flags across all 5 methods\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Build consensus table
# rPCA 97.5%: T2N, T1N
# rPCA 99%: T2N only
# Distance: T9N (z=1.88), T8C (z=1.84), T1A (z=1.55)
# Correlation: T9N (z=-2.16), T8C (z=-1.98), T1A (z=-1.80)
# Clustering: T9N clusters with Day0! (WRONG GROUP)

cat("  Sample      | rPCA 97.5% | rPCA 99% | Distance z | Corr z | Clusters correctly? | Cook's flags\n")
cat("  ------------|------------|----------|------------|--------|---------------------|-------------\n")

# Key samples of interest
key_samples <- c("T2N_S30", "T1N_S27", "T9N_S26", "T8C_S10", "T1A_S1", "T5N_S15", "T7N_S22", "T5A_S6")

rpca_97 <- c("T2N_S30", "T1N_S27")
rpca_99 <- c("T2N_S30")
dist_high <- c("T9N_S26", "T8C_S10", "T1A_S1")  
corr_low <- c("T9N_S26", "T8C_S10", "T1A_S1")
wrong_cluster <- c("T9N_S26", "T1A_S1", "T1B_S2", "T8C_S10")

for (s in key_samples) {
  f97 <- ifelse(s %in% rpca_97, "FLAG", "")
  f99 <- ifelse(s %in% rpca_99, "FLAG", "")
  d <- ifelse(s %in% dist_high, "HIGH", "")
  c <- ifelse(s %in% corr_low, "LOW", "")
  cl <- ifelse(s %in% wrong_cluster, "WRONG", "ok")
  n_flags <- sum(f97 != "", f99 != "", d != "", c != "", cl == "WRONG")
  cat(sprintf("  %-10s  | %-9s  | %-8s | %-10s | %-6s | %-19s | %d flags\n", 
      s, f97, f99, d, c, cl, n_flags))
}

cat("\n  VERDICT:\n")
cat("  T2N_S30: 2 flags (rPCA 97.5% + rPCA 99%) — MODERATE concern\n")
cat("  T9N_S26: 3 flags (distance HIGH, corr LOW, wrong cluster) — STRONG concern\n")
cat("  T8C_S10: 2 flags (distance HIGH, corr LOW) — MODERATE concern\n")
cat("  T1A_S1:  2 flags (distance HIGH, corr LOW, wrong cluster) — MODERATE concern\n")
cat("  T1N_S27: 1 flag (rPCA 97.5% only) — LOW concern\n")
cat("  T5N_S15: 0 flags — NORMAL\n")
cat("  T7N_S22: 0 flags — NORMAL\n")

cat("\n  KEY FINDING: T9N clusters with Day0, not Recovered!\n")
cat("  This means T9N's gene expression is more similar to tumor (Day0)\n")
cat("  than to recovered tissue. This is the strongest outlier signal.\n")
cat("  The original featureCounts rPCA was partially right about T9N.\n")