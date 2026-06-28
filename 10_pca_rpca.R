#!/usr/bin/env Rscript
# ====================================================================
# CTVT RNA-seq: Full 36-sample PCA + rPCA (Tumor + Recovered/Normal)
# Uses featureCounts gene-level data + deseq2_env
# ====================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(rrcov)
  library(robustbase)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(matrixStats)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
OUTDIR <- file.path(BASE, "rpca_deseq2_results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════════\n")
cat("  CTVT Full 36-sample PCA + rPCA Analysis\n")
cat("═══════════════════════════════════════════\n")

# ── Load data ──────────────────────────────────────────────────────
cat("\n[1/4] Loading featureCounts...\n")
fc_file <- file.path(BASE, "featureCounts_v3/gene_counts_v3.txt")
raw <- read.table(fc_file, header=TRUE, row.names=1, comment.char="#", check.names=FALSE)
counts_raw <- raw[, -(1:5)]

# Clean colnames
colnames(counts_raw) <- sub("\\.bam$", "", basename(colnames(counts_raw)))
colnames(counts_raw) <- sub("_Aligned\\.sortedByCoord\\.out$", "", colnames(counts_raw))

cat(sprintf("  %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# ── Build phenotype ────────────────────────────────────────────────
cat("\n[2/4] Building phenotype...\n")
extract_meta <- function(sid) {
  parts <- strsplit(sid, "_")[[1]]
  dog_id <- substr(parts[1], 1, 2)
  tissue_code <- substr(parts[1], 3, 3)
  sample_type <- ifelse(tissue_code == "N", "Normal", "Tumor")
  tp <- c(A="Day0", B="Day2", C="Day6", N="Recovered")[tissue_code]
  data.frame(sample_id=sid, dog_id=dog_id, timepoint=tp, 
             tissue_type=sample_type, tissue_code=tissue_code,
             stringsAsFactors=FALSE)
}

pheno <- do.call(rbind, lapply(colnames(counts_raw), extract_meta))
rownames(pheno) <- pheno$sample_id
counts_raw <- counts_raw[, pheno$sample_id]
stopifnot(all(colnames(counts_raw) == pheno$sample_id))

cat(sprintf("  %d samples: %d Tumor, %d Normal\n", nrow(pheno),
            sum(pheno$tissue_type=="Tumor"), sum(pheno$tissue_type=="Normal")))

# ── Pre-filter ──────────────────────────────────────────────────────
cat("\n[3/4] Filtering + VST...\n")
keep <- rowSums(counts_raw >= 10) >= 6
counts_filt <- counts_raw[keep, ]
cat(sprintf("  Kept %d / %d genes\n", nrow(counts_filt), nrow(counts_raw)))

# DESeq2 VST on all samples (blind design)
pheno$timepoint <- factor(pheno$timepoint, levels=c("Day0","Day2","Day6","Recovered"))
pheno$dog_id <- factor(pheno$dog_id)
pheno$tissue_type <- factor(pheno$tissue_type, levels=c("Tumor","Normal"))

dds_all <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = pheno,
  design = ~ 1
)
dds_all <- estimateSizeFactors(dds_all)
vsd <- vst(dds_all, blind=TRUE)
log_expr <- assay(vsd)

cat(sprintf("  VST matrix: %d x %d\n", nrow(log_expr), ncol(log_expr)))

# ── PCA + rPCA ──────────────────────────────────────────────────────
cat("\n[4/4] Running PCA + rPCA...\n")

# Top 500 most variable genes
rv <- rowVars(log_expr)
select_500 <- head(order(rv, decreasing=TRUE), 500)

# ── Standard PCA ────────────────────────────────────────────────────
pca_data <- t(log_expr[select_500, ])
pca_std <- prcomp(pca_data, center=TRUE, scale.=FALSE)
var_exp <- pca_std$sdev^2 / sum(pca_std$sdev^2) * 100
cat(sprintf("  Standard PCA: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%\n",
            var_exp[1], var_exp[2], var_exp[3]))

pca_df <- as.data.frame(pca_std$x)
pca_df$sample <- rownames(pca_df)
pca_df$dog <- pheno$dog_id
pca_df$timepoint <- pheno$timepoint
pca_df$tissue <- pheno$tissue_type

# ── Robust PCA (PcaGrid) ────────────────────────────────────────────
rpca <- PcaGrid(pca_data, k=5, scale=FALSE,
                method="mad", crit.pca.distances=0.975)
r_var_exp <- rpca$eigenvalues / sum(rpca$eigenvalues) * 100
cat(sprintf("  Robust PCA: rPC1=%.1f%%, rPC2=%.1f%%, rPC3=%.1f%%\n",
            r_var_exp[1], r_var_exp[2], r_var_exp[3]))

# Outlier detection
od <- rpca$od; sd <- rpca$sd
cutoff_od <- rpca$cutoff.od; cutoff_sd <- rpca$cutoff.sd

outlier_type <- rep("Normal", nrow(pca_data))
outlier_type[od > cutoff_od] <- "Orthogonal outlier"
outlier_type[sd > cutoff_sd] <- "Score outlier"
outlier_type[od > cutoff_od & sd > cutoff_sd] <- "Both (bad leverage)"

outliers_found <- which(outlier_type != "Normal")
cat(sprintf("  Outliers: %d / 36 samples\n", length(outliers_found)))
for (i in outliers_found) {
  cat(sprintf("    %s (%s, %s): %s (OD=%.1f cutoff=%.1f, SD=%.1f cutoff=%.1f)\n",
              rownames(pca_data)[i], pheno$tissue_type[i], pheno$timepoint[i],
              outlier_type[i], od[i], cutoff_od, sd[i], cutoff_sd))
}

rpca_df <- as.data.frame(rpca$scores)
colnames(rpca_df) <- paste0("rPC", 1:ncol(rpca_df))
rpca_df$sample <- rownames(pca_data)
rpca_df$dog <- pheno$dog_id
rpca_df$timepoint <- pheno$timepoint
rpca_df$tissue <- pheno$tissue_type
rpca_df$outlier <- outlier_type
rpca_df$od <- od
rpca_df$sd <- sd

# ══════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════
tp_colors <- c(Day0="#ef4444", Day2="#f97316", Day6="#eab308", Recovered="#22c55e")
dog_colors <- c(T1="#6366f1", T2="#8b5cf6", T3="#a855f7", T4="#d946ef",
                T5="#ec4899", T6="#f43f5e", T7="#14b8a6", T8="#0ea5e9", T9="#84cc16")

# ── PCA: Timepoint color, tissue shape ──
p1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=timepoint, shape=tissue)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=dog), size=2.8, max.overlaps=30,
                  box.padding=0.5, point.padding=0.3, show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  scale_shape_manual(values=c(Tumor=16, Normal=17)) +
  labs(title="CTVT: Standard PCA — All 36 Samples (Tumor + Normal)",
       subtitle="FeatureCounts gene-level, 500 MVG, VST-normalized",
       x=sprintf("PC1 (%.1f%%)", var_exp[1]),
       y=sprintf("PC2 (%.1f%%)", var_exp[2]),
       color="Timepoint", shape="Tissue") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"))

ggsave(file.path(OUTDIR, "PCA_full36_timepoint.png"), p1, width=11, height=8, dpi=150)

# ── PCA: Dog color, tissue shape ──
p2 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=dog, shape=tissue)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=timepoint), size=2.5, max.overlaps=30,
                  box.padding=0.5, show.legend=FALSE) +
  scale_color_manual(values=dog_colors) +
  scale_shape_manual(values=c(Tumor=16, Normal=17)) +
  labs(title="CTVT: PCA by Dog — All 36 Samples",
       x=sprintf("PC1 (%.1f%%)", var_exp[1]),
       y=sprintf("PC2 (%.1f%%)", var_exp[2]),
       color="Dog", shape="Tissue") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5))

ggsave(file.path(OUTDIR, "PCA_full36_by_dog.png"), p2, width=11, height=8, dpi=150)

# ── PCA: PC3 vs PC4 ──
p3 <- ggplot(pca_df, aes(x=PC3, y=PC4, color=timepoint, shape=tissue)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=dog), size=2.8, max.overlaps=30,
                  box.padding=0.5, show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  scale_shape_manual(values=c(Tumor=16, Normal=17)) +
  labs(title="CTVT: PC3 vs PC4",
       x=sprintf("PC3 (%.1f%%)", var_exp[3]),
       y=sprintf("PC4 (%.1f%%)", var_exp[4]),
       color="Timepoint", shape="Tissue") +
  theme_minimal(base_size=13)

ggsave(file.path(OUTDIR, "PCA_full36_pc3pc4.png"), p3, width=11, height=8, dpi=150)

# ── rPCA: Outlier plot ──
p4 <- ggplot(rpca_df, aes(x=rPC1, y=rPC2, color=timepoint, shape=outlier)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=paste0(dog, " (", tissue, ")")), 
                  size=2.5, max.overlaps=30, box.padding=0.5,
                  show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  scale_shape_manual(values=c(Normal=16, "Orthogonal outlier"=8,
                               "Score outlier"=15, "Both (bad leverage)"=17)) +
  labs(title="CTVT: Robust PCA (PcaGrid) — All 36 Samples",
       subtitle=sprintf("rPC1=%.1f%%, rPC2=%.1f%%  •  ★ = Outlier  •  (T)=Tumor (N)=Normal",
                        r_var_exp[1], r_var_exp[2]),
       x=sprintf("rPC1 (%.1f%%)", r_var_exp[1]),
       y=sprintf("rPC2 (%.1f%%)", r_var_exp[2]),
       color="Timepoint", shape="Outlier") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"))

ggsave(file.path(OUTDIR, "rPCA_full36_outliers.png"), p4, width=13, height=8, dpi=150)

# ── rPCA: Tissue separation ──
p5 <- ggplot(rpca_df, aes(x=rPC1, y=rPC2, color=tissue, shape=outlier)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=paste0(dog, "_", timepoint)),
                  size=2.5, max.overlaps=30, box.padding=0.4,
                  show.legend=FALSE) +
  scale_color_manual(values=c(Tumor="#ef4444", Normal="#22c55e")) +
  scale_shape_manual(values=c(Normal=16, "Orthogonal outlier"=8,
                               "Score outlier"=15, "Both (bad leverage)"=17)) +
  labs(title="CTVT: rPCA — Tumor vs Normal Separation",
       x=sprintf("rPC1 (%.1f%%)", r_var_exp[1]),
       y=sprintf("rPC2 (%.1f%%)", r_var_exp[2]),
       color="Tissue", shape="Outlier") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5))

ggsave(file.path(OUTDIR, "rPCA_full36_tissue.png"), p5, width=12, height=8, dpi=150)

# ── Diagnostic ──
png(file.path(OUTDIR, "rPCA_full36_diagnostic.png"), 
    width=10, height=8, units="in", res=150)
plot(rpca, main="PcaGrid Diagnostic: All 36 Samples")
dev.off()

# ── Scree plot ──
png(file.path(OUTDIR, "scree_full36.png"), width=12, height=6, units="in", res=150)
par(mfrow=c(1,2))
barplot(var_exp[1:10], names.arg=paste0("PC",1:10), col="#3b82f6",
        main="Standard PCA", ylab="% Variance")
text(seq_along(var_exp[1:10]), var_exp[1:10]+0.5,
     labels=sprintf("%.1f%%", var_exp[1:10]), cex=0.8)
barplot(r_var_exp[1:10], names.arg=paste0("rPC",1:10), col="#ef4444",
        main="Robust PCA (PcaGrid)", ylab="% Variance")
text(seq_along(r_var_exp[1:10]), r_var_exp[1:10]+0.5,
     labels=sprintf("%.1f%%", r_var_exp[1:10]), cex=0.8)
dev.off()

# ── Sample distance heatmap ──
sampleDists <- dist(t(log_expr))
annotation_col <- data.frame(
  Timepoint = pheno$timepoint,
  Dog = pheno$dog_id,
  Tissue = pheno$tissue_type,
  row.names = colnames(log_expr)
)
ann_colors <- list(
  Timepoint = tp_colors,
  Tissue = c(Tumor="#ef4444", Normal="#22c55e"),
  Dog = dog_colors
)
png(file.path(OUTDIR, "sample_distance_full36.png"),
    width=14, height=12, units="in", res=150)
pheatmap(as.matrix(sampleDists),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_distance_cols = sampleDists,
         clustering_method = "ward.D2",
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 6, fontsize_col = 6,
         color = colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(100),
         main = "CTVT: Sample Distance — All 36 Samples")
dev.off()

# ── Save results ──
write.csv(rpca_df, file.path(OUTDIR, "rpca_full36_outliers.csv"), row.names=FALSE)
write.csv(pca_df, file.path(OUTDIR, "pca_full36_scores.csv"), row.names=FALSE)

# ══════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════\n")
cat("  ✅ Full 36-sample analysis complete\n")
cat("═══════════════════════════════════════════\n")
cat(sprintf("  Results: %s/\n", OUTDIR))

cat("\n── TELEGRAM SUMMARY ──\n")
cat("📊 CTVT Full 36-Sample PCA + rPCA\n")
cat(sprintf("  Data: featureCounts gene-level, %d genes after filter\n", nrow(log_expr)))
cat(sprintf("  Samples: 27 Tumor (Day0/2/6) + 9 Normal (Recovered)\n"))
cat(sprintf("  Standard PCA: PC1=%.1f%%, PC2=%.1f%%\n", var_exp[1], var_exp[2]))
cat(sprintf("  Robust PCA: rPC1=%.1f%%, rPC2=%.1f%%\n", r_var_exp[1], r_var_exp[2]))
cat(sprintf("  rPCA outliers: %d / 36\n", length(outliers_found)))
if (length(outliers_found) > 0) {
  for (i in outliers_found) {
    cat(sprintf("    • %s (%s, %s): %s\n",
                rownames(pca_data)[i], pheno$tissue_type[i],
                pheno$timepoint[i], outlier_type[i]))
  }
}
cat(sprintf("\n  Recommendation: Do NOT auto-exclude outliers.\n"))
cat(sprintf("  Run DESeq2 with & without them, compare results.\n"))
cat(sprintf("  If DEG lists are stable → keep all samples.\n"))
cat(sprintf("  If outliers drive results → use rPCA weights or remove.\n"))
