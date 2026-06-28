#!/usr/bin/env Rscript
# CTVT: rPCA — Salmon + biotype-filtered + ~ dog_id + timepoint + removeBatchEffect(dog_id)
# This properly removes dog-to-dog variation before outlier detection
suppressPackageStartupMessages({
  library(DESeq2)
  library(rrcov)
  library(limma)
  library(tximport)
  library(matrixStats)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"

# Load biotype tx2gene (protein-coding + known lncRNA)
tx2gene_pc <- read.csv(file.path(BASE, "tx2gene_protein_coding.csv"), stringsAsFactors = FALSE)
tx2gene_lnc <- read.csv(file.path(BASE, "tx2gene_lncrna_known.csv"), stringsAsFactors = FALSE)
tx2gene_biotype <- rbind(tx2gene_pc, tx2gene_lnc)
cat(sprintf("Biotype tx2gene: %d transcripts\n", nrow(tx2gene_biotype)))

# Load ALL 36 samples
all_samples <- basename(list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE))
cat(sprintf("Total samples: %d\n", length(all_samples)))

sample_meta <- data.frame(sampleID = all_samples, stringsAsFactors = FALSE)
sample_meta$timepoint <- factor(sapply(sample_meta$sampleID, function(s) {
  c(A="Day0",B="Day2",C="Day6",N="Recovered")[substr(s,3,3)]
}), levels = c("Day0", "Day2", "Day6", "Recovered"))
sample_meta$dog_id <- factor(sapply(sample_meta$sampleID, function(s) {
  substr(strsplit(s, "_")[[1]][1], 1, 2)
}))

files <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
names(files) <- sample_meta$sampleID

# tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene_biotype, 
                  countsFromAbundance = "lengthScaledTPM")

# Filter
keep <- rowSums(txi$counts >= 10) >= 4
counts_raw <- round(txi$counts[keep, ])
cat(sprintf("After filter: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# DESeq2 with ~ dog_id + timepoint
dds <- DESeqDataSetFromMatrix(counts_raw, sample_meta, design = ~ dog_id + timepoint)
dds <- DESeq(dds, quiet = TRUE)

# VST (blind=FALSE: uses design for dispersion estimation)
vsd <- vst(dds, blind = FALSE)
log_expr <- assay(vsd)
cat(sprintf("VST matrix: %d genes x %d samples\n", nrow(log_expr), ncol(log_expr)))

# ══════════════════════════════════════════════════════════════════
# KEY STEP: Remove dog batch effect before PCA
# This ensures outlier detection sees residual variation AFTER
# accounting for the paired longitudinal design
# ══════════════════════════════════════════════════════════════════
log_expr_corrected <- removeBatchEffect(log_expr, batch = sample_meta$dog_id)
cat("\n── Removed dog_id batch effect via limma::removeBatchEffect ──\n")

# Top 500 variable genes from CORRECTED matrix
rv <- rowVars(log_expr_corrected)
select_500 <- head(order(rv, decreasing = TRUE), 500)
pca_data <- t(log_expr_corrected[select_500, ])

# Standard PCA on corrected data
pca_std <- prcomp(pca_data, center = TRUE, scale. = FALSE)
var_exp <- pca_std$sdev^2 / sum(pca_std$sdev^2) * 100
cat(sprintf("\n── Standard PCA (after removing dog batch) ──\n"))
cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n", var_exp[1], var_exp[2], var_exp[3]))

# Robust PCA on corrected data
cat("\n── Robust PCA (PcaGrid on batch-corrected data) ──\n")
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

results <- results[order(-pmax(results$OD_ratio, results$SD_ratio)), ]

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  rPCA — Salmon + biotype + ~dog_id+timepoint + removeBatchEffect(dog)\n")
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
cat("\n── Per-timepoint ──\n")
for (tp in c("Day0", "Day2", "Day6", "Recovered")) {
  idx <- results$Timepoint == tp
  n_flagged <- sum(results$Flag[idx] != "Normal")
  cat(sprintf("  %s: mean_OD=%.2fx, mean_SD=%.2fx, flagged=%d/%d\n", tp, 
      mean(results$OD_ratio[idx]), mean(results$SD_ratio[idx]),
      n_flagged, sum(idx)))
}

# Compare: also run WITHOUT removeBatchEffect for comparison
cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  COMPARISON: Without removeBatchEffect (raw VST)\n")
cat("═══════════════════════════════════════════════════════════════\n")
rv_raw <- rowVars(log_expr)
select_500_raw <- head(order(rv_raw, decreasing = TRUE), 500)
pca_data_raw <- t(log_expr[select_500_raw, ])

rpca_raw <- PcaGrid(pca_data_raw, k = 5, scale = FALSE, method = "mad", crit.pca.distances = 0.975)

od_raw <- rpca_raw$od
sd_raw <- rpca_raw$sd
cutoff_od_raw <- rpca_raw$cutoff.od
cutoff_sd_raw <- rpca_raw$cutoff.sd

outlier_raw <- rep("Normal", nrow(pca_data_raw))
outlier_raw[od_raw > cutoff_od_raw] <- "Orthogonal outlier"
outlier_raw[sd_raw > cutoff_sd_raw] <- "Score outlier"
outlier_raw[od_raw > cutoff_od_raw & sd_raw > cutoff_sd_raw] <- "Both (bad leverage)"

results_raw <- data.frame(
  Sample = rownames(pca_data_raw),
  OD_ratio_raw = round(od_raw / cutoff_od_raw, 2),
  SD_ratio_raw = round(sd_raw / cutoff_sd_raw, 2),
  Flag_raw = outlier_raw,
  stringsAsFactors = FALSE
)

# Merge
results_merged <- merge(results, results_raw, by = "Sample")
results_merged$changed <- results_merged$Flag != results_merged$Flag_raw
cat(sprintf("\n  Samples that CHANGED status after batch correction: %d\n", sum(results_merged$changed)))
cat(sprintf("  Before batch correction: %d flagged\n", sum(results_raw$Flag_raw != "Normal")))
cat(sprintf("  After batch correction:  %d flagged\n", sum(results$Flag != "Normal")))

for (s in original_outliers) {
  idx_c <- which(results$Sample == s)
  idx_r <- which(results_raw$Sample == s)
  if (length(idx_c) > 0 && length(idx_r) > 0) {
    cat(sprintf("\n  %s: CORRECTED=%.2fx/%.2fx(%s)  RAW=%.2fx/%.2fx(%s)", 
        s, results$OD_ratio[idx_c], results$SD_ratio[idx_c], results$Flag[idx_c],
        results_raw$OD_ratio_raw[idx_r], results_raw$SD_ratio_raw[idx_r], results_raw$Flag_raw[idx_r]))
  }
}

write.csv(results, file.path(BASE, "rpca_report_final_corrected.csv"), row.names = FALSE)
cat("\n\nSaved to rpca_report_final_corrected.csv\n")