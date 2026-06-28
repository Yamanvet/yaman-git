#!/usr/bin/env Rscript
# CTVT: rPCA with Salmon tximport + ~ dog_id + timepoint
# Using biotype-filtered genes (protein-coding + known lncRNA) to match our DESeq2 analysis
suppressPackageStartupMessages({
  library(DESeq2)
  library(rrcov)
  library(tximport)
  library(matrixStats)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"

# Load ALL 36 samples
all_samples <- basename(list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE))
cat(sprintf("Total samples: %d\n", length(all_samples)))

# Build metadata
sample_meta <- data.frame(sampleID = all_samples, stringsAsFactors = FALSE)
sample_meta$timepoint <- factor(sapply(sample_meta$sampleID, function(s) {
  c(A="Day0",B="Day2",C="Day6",N="Recovered")[substr(s,3,3)]
}), levels = c("Day0", "Day2", "Day6", "Recovered"))
sample_meta$dog_id <- factor(sapply(sample_meta$sampleID, function(s) {
  substr(strsplit(s, "_")[[1]][1], 1, 2)
}))

# Load tximport
tx2gene_file <- file.path(BASE, "salmon_index/tx2gene.csv")
tx2gene <- read.csv(tx2gene_file, stringsAsFactors = FALSE)
files <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
names(files) <- sample_meta$sampleID

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                  countsFromAbundance = "lengthScaledTPM")

# Filter: >=10 counts in >=4 samples
keep <- rowSums(txi$counts >= 10) >= 4
counts_raw <- round(txi$counts[keep, ])
cat(sprintf("After low-count filter: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# ── DESeq2 with ~ dog_id + timepoint ──
dds <- DESeqDataSetFromMatrix(counts_raw, sample_meta, design = ~ dog_id + timepoint)
cat(sprintf("DESeq2 design: ~ dog_id + timepoint\n"))
dds <- DESeq(dds, quiet = TRUE)

# VST (blind=FALSE: uses the design to remove dog effect)
cat("\n── VST with design ~ dog_id + timepoint ──\n")
vsd <- vst(dds, blind = FALSE)
log_expr <- assay(vsd)
cat(sprintf("  VST matrix: %d genes x %d samples\n", nrow(log_expr), ncol(log_expr)))

# ── Robust PCA on top 500 variable genes ──
rv <- rowVars(log_expr)
select_500 <- head(order(rv, decreasing = TRUE), 500)
pca_data <- t(log_expr[select_500, ])

# Standard PCA first
pca_std <- prcomp(pca_data, center = TRUE, scale. = FALSE)
var_exp <- pca_std$sdev^2 / sum(pca_std$sdev^2) * 100
cat(sprintf("\n── Standard PCA (with dog_id+timepoint VST) ──\n"))
cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n", var_exp[1], var_exp[2], var_exp[3]))

# Robust PCA
cat("\n── Robust PCA (PcaGrid, k=5, 97.5% cutoff) ──\n")
rpca <- PcaGrid(pca_data, k = 5, scale = FALSE, method = "mad", crit.pca.distances = 0.975)

od <- rpca$od
sd_val <- rpca$sd
cutoff_od <- rpca$cutoff.od
cutoff_sd <- rpca$cutoff.sd
cat(sprintf("  OD cutoff: %.2f | SD cutoff: %.2f\n", cutoff_od, cutoff_sd))

# Build report
outlier_flags <- rep("Normal", nrow(pca_data))
outlier_flags[od > cutoff_od] <- "Orthogonal outlier"
outlier_flags[sd_val > cutoff_sd] <- "Score outlier"
outlier_flags[od > cutoff_od & sd_val > cutoff_sd] <- "Both (bad leverage)"

results <- data.frame(
  Sample = rownames(pca_data),
  Timepoint = sample_meta$timepoint[match(rownames(pca_data), sample_meta$sampleID)],
  Dog = sample_meta$dog_id[match(rownames(pca_data), sample_meta$sampleID)],
  OD = round(od, 2),
  OD_cutoff = round(cutoff_od, 2),
  OD_ratio = round(od / cutoff_od, 2),
  SD = round(sd_val, 2),
  SD_cutoff = round(cutoff_sd, 2),
  SD_ratio = round(sd_val / cutoff_sd, 2),
  Flag = outlier_flags,
  stringsAsFactors = FALSE
)

# Sort by worst
results <- results[order(-pmax(results$OD_ratio, results$SD_ratio)), ]

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  rPCA REPORT — Salmon + ~ dog_id + timepoint (blind=FALSE)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (i in 1:nrow(results)) {
  r <- results[i, ]
  marker <- ifelse(r$Flag != "Normal", "!! ", "   ")
  cat(sprintf("%s %-8s Dog=%s Time=%-10s  OD=%5.1f (%.2fx)  SD=%4.1f (%.2fx)  %s\n",
      marker, r$Sample, r$Dog, r$Timepoint,
      r$OD, r$OD_ratio, r$SD, r$SD_ratio, r$Flag))
}

flagged <- sum(outlier_flags != "Normal")
cat(sprintf("\n\n── Summary ──\n"))
cat(sprintf("  Flagged: %d / %d samples\n", flagged, length(all_samples)))

# Check original outliers
original_outliers <- c("T1B_S2", "T8C_S10", "T9N_S26")
cat("\n── Original outlier positions ──\n")
for (s in original_outliers) {
  idx <- which(results$Sample == s)
  if (length(idx) > 0) {
    r <- results[idx, ]
    cat(sprintf("  %s: OD=%.1f (%.2fx), SD=%.1f (%.2fx), Rank=%d/%d, Flag=%s\n",
        s, r$OD, r$OD_ratio, r$SD, r$SD_ratio, idx, nrow(results), r$Flag))
  }
}

# Per-timepoint
cat("\n── Per-timepoint mean OD/SD ratios ──\n")
for (tp in c("Day0", "Day2", "Day6", "Recovered")) {
  idx <- results$Timepoint == tp
  n_flagged <- sum(results$Flag[idx] != "Normal")
  cat(sprintf("  %s: mean_OD=%.2fx, mean_SD=%.2fx, flagged=%d/%d\n", tp, 
      mean(results$OD_ratio[idx]), mean(results$SD_ratio[idx]),
      n_flagged, sum(idx)))
}

write.csv(results, file.path(BASE, "rpca_report_salmon_dog_model.csv"), row.names = FALSE)
cat("\nSaved to rpca_report_salmon_dog_model.csv\n")