#!/usr/bin/env Rscript
# ====================================================================
# CTVT RNA-seq: Robust PCA + DESeq2 Differential Expression
# Uses deseq2_env conda environment (R 4.3.3, DESeq2 1.42.0)
# Memory: R_MAX_VSIZE=250G (server has 1TB, leave headroom)
# ====================================================================

# ── Setup ──────────────────────────────────────────────────────────
Sys.setenv(R_MAX_VSIZE = "250e9")  # 250 GB max per session
options(future.globals.maxSize = 10 * 1024^3)

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
cat("  CTVT rPCA + DESeq2 Analysis Pipeline\n")
cat("═══════════════════════════════════════════\n")

# ── Load featureCounts gene-level data ─────────────────────────────
cat("\n[1/6] Loading featureCounts gene counts...\n")
fc_file <- file.path(BASE, "featureCounts_v3/gene_counts_v3.txt")
raw <- read.table(fc_file, header=TRUE, row.names=1, comment.char="#", check.names=FALSE)

# Columns 1-5 are annotation (Chr, Start, End, Strand, Length), rest are samples
gene_annot <- raw[, 1:5]
counts_raw <- raw[, -(1:5)]

# Clean sample names
colnames(counts_raw) <- sub("\\.bam$", "", basename(colnames(counts_raw)))
colnames(counts_raw) <- sub("_Aligned\\.sortedByCoord\\.out$", "", colnames(counts_raw))
stopifnot(all(grepl("^T[0-9][A-Z]_S[0-9]+$", colnames(counts_raw))))

cat(sprintf("  Loaded: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# ── Build phenotype ────────────────────────────────────────────────
cat("\n[2/6] Building phenotype table...\n")
extract_meta <- function(sid) {
  parts <- strsplit(sid, "_")[[1]]
  # e.g. T1A_S1 -> dog=T1, tissue=A (Tumor), replicate=S1
  dog_id <- substr(parts[1], 1, 2)  # T1, T2, ... T9
  tissue_code <- substr(parts[1], 3, 3)
  sample_type <- ifelse(tissue_code == "N", "Normal", "Tumor")
  # Timepoint from tissue_code: A=Day0, B=Day2, C=Day6, N=Recovered
  timepoint_map <- c(A="Day0", B="Day2", C="Day6", N="Recovered")
  timepoint <- timepoint_map[tissue_code]
  data.frame(sample_id=sid, dog_id=dog_id, timepoint=timepoint, 
             tissue_type=sample_type, tissue_code=tissue_code,
             stringsAsFactors=FALSE)
}

pheno <- do.call(rbind, lapply(colnames(counts_raw), extract_meta))
rownames(pheno) <- pheno$sample_id

# Reorder columns to match phenotype
counts_raw <- counts_raw[, pheno$sample_id]
stopifnot(all(colnames(counts_raw) == pheno$sample_id))

cat(sprintf("  Samples: %d total (%d Tumor, %d Normal)\n",
            nrow(pheno), sum(pheno$tissue_type=="Tumor"), sum(pheno$tissue_type=="Normal")))
cat(sprintf("  Dogs: %d\n", length(unique(pheno$dog_id))))
cat("  Timepoints: Day0, Day2, Day6, Recovered\n")

# ── Pre-filter low-count genes ─────────────────────────────────────
cat("\n[3/6] Filtering low-count genes...\n")
# Keep genes with >= 10 counts in >= 4 samples (flexible for small groups)
keep <- rowSums(counts_raw >= 10) >= 4
counts_filt <- counts_raw[keep, ]
gene_annot_filt <- gene_annot[keep, ]
cat(sprintf("  Kept %d / %d genes after filtering\n", nrow(counts_filt), nrow(counts_raw)))

# ── DESeq2: All tumor samples (timepoint design) ───────────────────
cat("\n[4/6] Running DESeq2 on tumor-only samples...\n")

# Tumor only for time course
tumor_idx <- which(pheno$tissue_type == "Tumor")
tumor_counts <- counts_filt[, tumor_idx, drop=FALSE]
tumor_pheno <- pheno[tumor_idx, ]
tumor_pheno$timepoint <- factor(tumor_pheno$timepoint, levels=c("Day0","Day2","Day6","Recovered"))
tumor_pheno$dog_id <- factor(tumor_pheno$dog_id)

dds <- DESeqDataSetFromMatrix(
  countData = tumor_counts,
  colData = tumor_pheno,
  design = ~ dog_id + timepoint
)

# NOTE: Already pre-filtered at gene level above.
# Do NOT re-filter inside DESeq2 — some timepoints (Recovered=Normal) have very
# different expression and would get dropped, breaking the model matrix.
# Just use the pre-filtered counts directly.
dds <- dds[rowSums(counts(dds)) > 0, ]  # minimal: drop truly dead genes

# Run DESeq2 (this is the heavy step)
cat("  Estimating size factors...\n")
dds <- estimateSizeFactors(dds)
cat("  Estimating dispersions...\n")
dds <- estimateDispersions(dds, quiet=FALSE)
cat("  Fitting nbinomWaldTest...\n")
dds <- nbinomWaldTest(dds, maxit=5000, useOptim=TRUE)

# ── DEG extraction: All contrasts vs Day0 ──────────────────────────
cat("\n[5/6] Extracting DEG results...\n")

# Note: Recovered=N samples are normal tissue, not tumor timepoint.
# Tumor-only analysis has only Day0, Day2, Day6.
# Get exact coefficient names from model
rn <- resultsNames(dds)
cat("  Available coefficients:\n")
for (cc in grep("timepoint", rn, value=TRUE)) cat(sprintf("    %s\n", cc))

deseq_results <- list()
contrast_map <- list(
  "Day2_vs_Day0"      = grep("timepoint.*Day2", rn, value=TRUE)[1],
  "Day6_vs_Day0"      = grep("timepoint.*Day6", rn, value=TRUE)[1],
  "Day6_vs_Day2"      = grep("timepoint.*Day6", rn, value=TRUE)[1]
)

for (cname in names(contrast_map)) {
  coef_name <- contrast_map[[cname]]
  if (is.na(coef_name) || !coef_name %in% rn) {
    cat(sprintf("  SKIPPING %s: coefficient not found\n", cname))
    next
  }
  cat(sprintf("  %s [%s]...\n", cname, coef_name))
  
  res <- results(dds, name=coef_name, alpha=0.05, 
                 lfcThreshold=0, independentFiltering=TRUE)
  res <- res[order(res$pvalue), ]
  
  # For Day6_vs_Day2, we need the difference of coefficients
  if (cname == "Day6_vs_Day2") {
    # LFC = Day6 - Day2; use contrast with Day2 as reference
    res <- results(dds, contrast=c("timepoint", "Day6", "Day2"), alpha=0.05)
  }
  
  res <- res[order(res$pvalue), ]
  
  # Add gene annotation
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # Significance
  res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "Significant", "NS")
  res_df$direction <- ifelse(res_df$significant == "Significant",
                              ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "NS")
  
  outfile <- file.path(OUTDIR, paste0("DEG_", cname, ".csv"))
  write.csv(res_df[order(res_df$pvalue), ], outfile, row.names=FALSE)
  cat(sprintf("    Saved: %s\n", basename(outfile)))
  
  sig_count <- sum(res_df$significant == "Significant", na.rm=TRUE)
  up_count <- sum(res_df$direction == "Up", na.rm=TRUE)
  down_count <- sum(res_df$direction == "Down", na.rm=TRUE)
  cat(sprintf("    %d DEGs (padj<0.05): %d up, %d down\n", sig_count, up_count, down_count))
  
  deseq_results[[cname]] <- res_df
}

# Also: Day6 vs Day2 contrast (already handled above in loop)
# No Recovered timepoint in tumor-only — Recovered/N are normal tissue samples
sig_r6 <- sum(deseq_results[["Day6_vs_Day2"]]$significant == "Significant", na.rm=TRUE)

# ── VST transformation for PCA ─────────────────────────────────────
cat("\n[6/6] VST transformation + PCA / rPCA...\n")

vsd <- vst(dds, blind=TRUE)
log_expr <- assay(vsd)
cat(sprintf("  VST matrix: %d genes x %d samples\n", nrow(log_expr), ncol(log_expr)))

# Remove batch effect? Check if dog explains a lot
# We'll do PCA on the VST residuals after removing dog effect
# For now, also do raw VST PCA for comparison

# ══════════════════════════════════════════════════════════════════
# STANDARD PCA (500 most variable genes)
# ══════════════════════════════════════════════════════════════════
cat("\n── Standard PCA ──\n")
rv <- rowVars(log_expr)
select_500 <- head(order(rv, decreasing=TRUE), 500)

pca_std <- prcomp(t(log_expr[select_500, ]), center=TRUE, scale.=FALSE)
var_exp <- pca_std$sdev^2 / sum(pca_std$sdev^2) * 100
cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n", var_exp[1], var_exp[2], var_exp[3]))

# PCA dataframe
pca_df <- as.data.frame(pca_std$x)
pca_df$sample <- rownames(pca_df)
pca_df$dog <- colData(dds)$dog_id
pca_df$timepoint <- colData(dds)$timepoint

# ══════════════════════════════════════════════════════════════════
# ROBUST PCA (PcaGrid - optimal for outlier detection)
# ══════════════════════════════════════════════════════════════════
cat("\n── Robust PCA (PcaGrid) ──\n")
# Use top 500 most variable genes, same as standard PCA
pca_data <- t(log_expr[select_500, ])

# PcaGrid: robust PCA using projection pursuit + grid search
# k=5 components to capture structure
rpca <- PcaGrid(pca_data, k=5, scale=FALSE, 
                method="mad", crit.pca.distances=0.975)

cat(sprintf("  Robust PCA eigenvalues: %.1f %.1f %.1f %.1f %.1f\n",
            rpca$eigenvalues[1], rpca$eigenvalues[2], rpca$eigenvalues[3],
            rpca$eigenvalues[4], rpca$eigenvalues[5]))

r_var_exp <- rpca$eigenvalues / sum(rpca$eigenvalues) * 100

# Flag outliers: orthogonal distance > cutoff
od <- rpca$od
sd <- rpca$sd
cutoff_od <- rpca$cutoff.od
cutoff_sd <- rpca$cutoff.sd

outlier_flags <- rep("Normal", nrow(pca_data))
outlier_flags[od > cutoff_od] <- "Orthogonal outlier"
outlier_flags[sd > cutoff_sd] <- "Score outlier" 
outlier_flags[od > cutoff_od & sd > cutoff_sd] <- "Both (bad leverage)"

cat(sprintf("  Outliers detected: %d samples flagged\n", sum(outlier_flags != "Normal")))
for (i in which(outlier_flags != "Normal")) {
  cat(sprintf("    %s: %s (OD=%.1f, SD=%.1f)\n", 
              rownames(pca_data)[i], outlier_flags[i], od[i], sd[i]))
}

# rPCA scores dataframe
rpca_df <- as.data.frame(rpca$scores)
colnames(rpca_df) <- paste0("rPC", 1:ncol(rpca_df))
rpca_df$sample <- rownames(pca_data)
rpca_df$dog <- colData(dds)$dog_id
rpca_df$timepoint <- colData(dds)$timepoint
rpca_df$outlier <- outlier_flags
rpca_df$od <- od
rpca_df$sd <- sd

# ══════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════

# Timepoint colors
tp_colors <- c(Day0="#ef4444", Day2="#f97316", Day6="#eab308", Recovered="#22c55e")
dog_colors <- c(T1="#6366f1", T2="#8b5cf6", T3="#a855f7", T4="#d946ef",
                T5="#ec4899", T6="#f43f5e", T7="#14b8a6", T8="#0ea5e9", T9="#84cc16")

# ── Standard PCA by timepoint ──
p1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=timepoint)) +
  geom_point(size=4.5, alpha=0.85) +
  geom_text_repel(aes(label=dog), size=3.2, max.overlaps=20, 
                  box.padding=0.5, point.padding=0.3, show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  stat_ellipse(aes(group=timepoint), type="norm", level=0.95, 
               linewidth=0.8, alpha=0.4, show.legend=FALSE) +
  labs(title="CTVT Tumor RNA-seq: Standard PCA (500 MVG)",
       subtitle="All timepoints - Tumor Only",
       x=sprintf("PC1 (%.1f%%)", var_exp[1]),
       y=sprintf("PC2 (%.1f%%)", var_exp[2]),
       color="Timepoint") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"),
        panel.grid.minor=element_blank())

ggsave(file.path(OUTDIR, "PCA_standard_timepoint.png"), p1, width=10, height=8, dpi=150)

# ── Standard PCA by dog ──
p2 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=dog)) +
  geom_point(size=4.5, alpha=0.85) +
  geom_text_repel(aes(label=timepoint), size=3, max.overlaps=20,
                  box.padding=0.5, show.legend=FALSE) +
  scale_color_manual(values=dog_colors) +
  labs(title="CTVT Tumor RNA-seq: PCA by Dog",
       subtitle="Colors = Dog, Labels = Timepoint",
       x=sprintf("PC1 (%.1f%%)", var_exp[1]),
       y=sprintf("PC2 (%.1f%%)", var_exp[2]),
       color="Dog") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"))

ggsave(file.path(OUTDIR, "PCA_standard_by_dog.png"), p2, width=10, height=8, dpi=150)

# ── rPCA score plot ──
p3 <- ggplot(rpca_df, aes(x=rPC1, y=rPC2, color=timepoint, shape=outlier)) +
  geom_point(size=4.5, alpha=0.85) +
  geom_text_repel(aes(label=dog), size=3.2, max.overlaps=20,
                  box.padding=0.5, show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  scale_shape_manual(values=c(Normal=16, "Orthogonal outlier"=8, 
                               "Score outlier"=15, "Both (bad leverage)"=17)) +
  labs(title="CTVT: Robust PCA (PcaGrid) — Outlier Detection",
       subtitle=sprintf("rPC1 (%.1f%%) vs rPC2 (%.1f%%) • ★ = outlier", 
                        r_var_exp[1], r_var_exp[2]),
       x=sprintf("rPC1 (%.1f%%)", r_var_exp[1]),
       y=sprintf("rPC2 (%.1f%%)", r_var_exp[2]),
       color="Timepoint", shape="Outlier Status") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"))

ggsave(file.path(OUTDIR, "rPCA_outliers.png"), p3, width=11, height=8, dpi=150)

# ── rPCA by dog ──
p4 <- ggplot(rpca_df, aes(x=rPC1, y=rPC2, color=dog)) +
  geom_point(aes(shape=outlier), size=4.5, alpha=0.85) +
  geom_text_repel(aes(label=timepoint), size=3, max.overlaps=20,
                  box.padding=0.5, show.legend=FALSE) +
  scale_color_manual(values=dog_colors) +
  scale_shape_manual(values=c(Normal=16, "Orthogonal outlier"=8,
                               "Score outlier"=15, "Both (bad leverage)"=17)) +
  labs(title="CTVT: Robust PCA by Dog",
       subtitle="★ = outlier detected by PcaGrid",
       x=sprintf("rPC1 (%.1f%%)", r_var_exp[1]),
       y=sprintf("rPC2 (%.1f%%)", r_var_exp[2]),
       color="Dog", shape="Outlier") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"))

ggsave(file.path(OUTDIR, "rPCA_by_dog.png"), p4, width=10, height=8, dpi=150)

# ── Distance diagnostic plot ──
png(file.path(OUTDIR, "rPCA_diagnostic.png"), width=10, height=8, units="in", res=150)
plot(rpca, main="PcaGrid: Outlier Diagnostic Plot")
dev.off()

# ── Scree plot comparison ──
png(file.path(OUTDIR, "scree_plot_comparison.png"), width=12, height=6, units="in", res=150)
par(mfrow=c(1,2))
# Standard PCA scree
barplot(var_exp[1:10], names.arg=paste0("PC",1:10), col="#3b82f6",
        main="Standard PCA: Variance Explained", ylab="% Variance", 
        ylim=c(0, max(var_exp[1:10])*1.15))
text(x=seq_along(var_exp[1:10]), y=var_exp[1:10]+0.5, 
     labels=sprintf("%.1f%%", var_exp[1:10]), cex=0.8)

# rPCA scree
barplot(r_var_exp[1:min(10,length(r_var_exp))], 
        names.arg=paste0("rPC",1:min(10,length(r_var_exp))), col="#ef4444",
        main="Robust PCA (PcaGrid): Variance Explained", ylab="% Variance",
        ylim=c(0, max(r_var_exp[1:min(10,length(r_var_exp))])*1.15))
text(x=seq_len(min(10,length(r_var_exp))), y=r_var_exp[1:min(10,length(r_var_exp))]+0.5,
     labels=sprintf("%.1f%%", r_var_exp[1:min(10,length(r_var_exp))]), cex=0.8)
dev.off()

# ══════════════════════════════════════════════════════════════════
# HEATMAP: Top DEGs across contrasts
# ══════════════════════════════════════════════════════════════════
cat("\n── Generating heatmaps ──\n")

# Combine top DEGs from all contrasts
top_genes_list <- list()
for (cname in names(deseq_results)) {
  res <- deseq_results[[cname]]
  sig_genes <- rownames(res)[which(res$significant == "Significant")]
  top_genes_list[[cname]] <- head(sig_genes, 100)
}

# Union of top 100 from each contrast, cap at 200
all_top <- unique(unlist(top_genes_list))
all_top <- head(all_top, 200)

cat(sprintf("  Heatmap genes: %d (union of top 100 from each contrast)\n", length(all_top)))

# Subset expression
heatmap_data <- log_expr[all_top[all_top %in% rownames(log_expr)], , drop=FALSE]

# Scale rows for visualization
heatmap_scaled <- t(scale(t(heatmap_data)))

# Clamp extremes for color scaling
heatmap_scaled[heatmap_scaled > 3] <- 3
heatmap_scaled[heatmap_scaled < -3] <- -3

# Annotation
annotation_col <- data.frame(
  Timepoint = colData(dds)$timepoint,
  Dog = colData(dds)$dog_id,
  row.names = colnames(heatmap_scaled)
)

ann_colors <- list(
  Timepoint = tp_colors,
  Dog = dog_colors
)

# Shorten gene names
gene_labels <- rownames(heatmap_scaled)
gene_labels <- substr(gene_labels, 1, 18)

png(file.path(OUTDIR, "heatmap_topDEGs.png"), 
    width=16, height=max(10, length(all_top)*0.25), units="in", res=150)
pheatmap(heatmap_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_row = gene_labels,
         show_rownames = if(nrow(heatmap_scaled) <= 100) TRUE else FALSE,
         show_colnames = TRUE,
         fontsize_row = 5,
         fontsize_col = 7,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "CTVT: Top DEGs (Union of All Contrasts)",
         border_color = NA)
dev.off()

# ══════════════════════════════════════════════════════════════════
# VOLCANO PLOTS per contrast
# ══════════════════════════════════════════════════════════════════
cat("\n── Volcano plots ──\n")

for (cname in names(deseq_results)) {
  res <- deseq_results[[cname]]
  res$negLog10P <- -log10(res$pvalue)
  
  # Top 10 genes to label
  top_genes <- head(rownames(res)[order(res$pvalue)], 10)
  res$label <- ifelse(rownames(res) %in% top_genes, 
                      substr(rownames(res), 1, 15), "")
  
  p <- ggplot(res, aes(x=log2FoldChange, y=negLog10P, color=direction)) +
    geom_point(size=1.5, alpha=0.6) +
    geom_text_repel(aes(label=label), size=3, max.overlaps=15,
                    box.padding=0.3, na.rm=TRUE) +
    scale_color_manual(values=c(Up="#ef4444", Down="#3b82f6", NS="grey70")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
    geom_vline(xintercept=c(-1,1), linetype="dotted", color="grey50") +
    labs(title=paste("CTVT:", gsub("_", " ", cname)),
         x="log2 Fold Change",
         y="-log10(p-value)",
         color="Direction") +
    theme_minimal(base_size=13) +
    theme(plot.title=element_text(face="bold", hjust=0.5))
  
  ggsave(file.path(OUTDIR, paste0("volcano_", cname, ".png")), 
         p, width=10, height=8, dpi=150)
  cat(sprintf("  %s saved\n", cname))
}

# ══════════════════════════════════════════════════════════════════
# SAMPLE DISTANCE HEATMAP
# ══════════════════════════════════════════════════════════════════
cat("\n── Sample distance heatmap ──\n")
sampleDists <- dist(t(log_expr))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(log_expr)
colnames(sampleDistMatrix) <- colnames(log_expr)

png(file.path(OUTDIR, "sample_distance_heatmap.png"), 
    width=12, height=10, units="in", res=150)
pheatmap(sampleDistMatrix,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         fontsize_row = 7,
         fontsize_col = 7,
         color = colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(100),
         main = "CTVT: Sample-to-Sample Euclidean Distance (VST)")
dev.off()

# ══════════════════════════════════════════════════════════════════
# DEG OVERLAP (Upset-style count)
# ══════════════════════════════════════════════════════════════════
cat("\n── DEG overlap summary ──\n")

deg_sets <- list()
for (cname in names(deseq_results)) {
  res <- deseq_results[[cname]]
  sig <- rownames(res)[which(res$significant == "Significant")]
  deg_sets[[cname]] <- sig
  cat(sprintf("  %s: %d DEGs\n", cname, length(sig)))
}

# Pairwise overlaps
for (i in 1:(length(deg_sets)-1)) {
  for (j in (i+1):length(deg_sets)) {
    overlap <- intersect(deg_sets[[i]], deg_sets[[j]])
    cat(sprintf("  %s ∩ %s: %d genes\n", names(deg_sets)[i], names(deg_sets)[j], length(overlap)))
  }
}

# 3-way overlap
three_way <- intersect(intersect(deg_sets[[1]], deg_sets[[2]]), deg_sets[[3]])
cat(sprintf("  All 3 contrasts: %d shared genes\n", length(three_way)))

# Save overlap info
overlap_df <- data.frame(
  contrast = names(deg_sets),
  n_DEGs = sapply(deg_sets, length),
  n_up = sapply(names(deg_sets), function(n) sum(deseq_results[[n]]$direction=="Up", na.rm=TRUE)),
  n_down = sapply(names(deg_sets), function(n) sum(deseq_results[[n]]$direction=="Down", na.rm=TRUE))
)
write.csv(overlap_df, file.path(OUTDIR, "DEG_summary.csv"), row.names=FALSE)

# ══════════════════════════════════════════════════════════════════
# SAVE FULL RESULTS
# ══════════════════════════════════════════════════════════════════
cat("\n── Saving full results ──\n")

# Save VST matrix for downstream
write.csv(as.data.frame(log_expr), file.path(OUTDIR, "vst_expression_matrix.csv"))
cat("  Saved: vst_expression_matrix.csv\n")

# Save rPCA results
write.csv(rpca_df, file.path(OUTDIR, "rpca_scores_outliers.csv"), row.names=FALSE)
cat("  Saved: rpca_scores_outliers.csv\n")

# Save normalized counts
norm_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(norm_counts), file.path(OUTDIR, "deseq2_normalized_counts.csv"))
cat("  Saved: deseq2_normalized_counts.csv\n")

# ══════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════
cat("\n═══════════════════════════════════════════\n")
cat("  ✅ Analysis Complete!\n")
cat("═══════════════════════════════════════════\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("\nFiles generated:\n")
for (f in sort(list.files(OUTDIR))) {
  cat(sprintf("  %s\n", f))
}

# Quick summary for Telegram
cat("\n── TELEGRAM SUMMARY ──\n")
cat(sprintf("📊 CTVT RNA-seq Analysis Complete\n\n"))
cat(sprintf("🟢 DESeq2: %d genes x %d tumor samples\n", nrow(counts(dds)), ncol(counts(dds))))
for (i in seq_along(deg_sets)) {
  cat(sprintf("  %s: %d DEGs (%d↑ %d↓)\n",
              names(deg_sets)[i], 
              overlap_df$n_DEGs[i], overlap_df$n_up[i], overlap_df$n_down[i]))
}
cat(sprintf("  Day6_vs_Day2: %d DEGs\n", sig_r6))
cat(sprintf("\n🔴 rPCA (PcaGrid): %d potential outliers detected\n", sum(outlier_flags != "Normal")))
if (sum(outlier_flags != "Normal") > 0) {
  for (i in which(outlier_flags != "Normal")) {
    cat(sprintf("  • %s: %s\n", rownames(pca_data)[i], outlier_flags[i]))
  }
}
cat(sprintf("\n📈 Standard PCA: PC1=%.1f%%, PC2=%.1f%%\n", var_exp[1], var_exp[2]))
cat(sprintf("📈 Robust PCA: rPC1=%.1f%%, rPC2=%.1f%%\n", r_var_exp[1], r_var_exp[2]))
cat(sprintf("\n📁 Results in: %s/\n", OUTDIR))
