#!/usr/bin/env Rscript
# CTVT: Full rPCA outlier report for all 36 samples
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

# Filter
keep <- rowSums(txi$counts >= 10) >= 4
counts_filt <- txi$counts[keep, ]
cat(sprintf("Filtered genes: %d\n", nrow(counts_filt)))

# Round to integers for DESeq2 (tximport lengthScaledTPM produces floats)
counts_filt <- round(counts_filt)

# VST transform
dds <- DESeqDataSetFromMatrix(counts_filt, sample_meta, ~ timepoint)
vsd <- vst(dds, blind = TRUE)
log_expr <- assay(vsd)

# Standard PCA
rv <- rowVars(log_expr)
select_500 <- head(order(rv, decreasing = TRUE), 500)
pca_data <- t(log_expr[select_500, ])

# Robust PCA (PcaGrid)
rpca <- PcaGrid(pca_data, k = 5, scale = FALSE, method = "mad", crit.pca.distances = 0.975)

od <- rpca$od
sd_val <- rpca$sd
cutoff_od <- rpca$cutoff.od
cutoff_sd <- rpca$cutoff.sd

cat(sprintf("\nCutoffs: OD=%.2f, SD=%.2f\n", cutoff_od, cutoff_sd))

# Build full report
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

# Sort by worst outlier
results <- results[order(-pmax(results$OD_ratio, results$SD_ratio)), ]

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  rPCA OUTLIER REPORT — ALL 36 SAMPLES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (i in 1:nrow(results)) {
  r <- results[i, ]
  marker <- ifelse(r$Flag != "Normal", "!! ", "   ")
  cat(sprintf("%s %s  Dog=%s Time=%-10s  OD=%.1f (%.2fx)  SD=%.1f (%.2fx)  %s\n",
      marker, r$Sample, r$Dog, r$Timepoint,
      r$OD, r$OD_ratio, r$SD, r$SD_ratio, r$Flag))
}

# Summary
cat(sprintf("\n\n-- Summary --\n"))
cat(sprintf("  OD cutoff: %.2f | SD cutoff: %.2f\n", cutoff_od, cutoff_sd))
cat(sprintf("  Flagged: %d / %d samples\n", sum(outlier_flags != "Normal"), length(all_samples)))

# Per-timepoint stats
cat("\n-- Per-timepoint mean OD/SD ratios --\n")
for (tp in c("Day0", "Day2", "Day6", "Recovered")) {
  idx <- results$Timepoint == tp
  cat(sprintf("  %s: mean_OD_ratio=%.2f, mean_SD_ratio=%.2f\n", tp, 
      mean(results$OD_ratio[idx]), mean(results$SD_ratio[idx])))
}

# Save
write.csv(results, file.path(BASE, "rpca_outlier_report_full.csv"), row.names = FALSE)
cat("\nSaved to rpca_outlier_report_full.csv\n")