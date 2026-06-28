#!/usr/bin/env Rscript
# ====================================================================
# CTVT: Isoform-level PCA + DTU (fast, practical)
# Uses transcript counts → known gene symbols via tx2gene_real
# PCA: Top 500 most variable transcripts (isoform-level)
# DTU: Per-gene proportion shift test (chi-squared, filtered)
# ====================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
OUTDIR <- file.path(BASE, "isoform_pca_dtu")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

cat("═══════════════════════════════════════════\n")
cat("  CTVT Isoform-Level PCA + DTU\n")
cat("═══════════════════════════════════════════\n")

# ── Load ───────────────────────────────────────────────────────────
cat("\n[1/5] Loading data...\n")

tx_counts <- read.csv(file.path(BASE, "ballgown_v3/transcript_count_matrix_prepDE.csv"),
                      row.names=1, check.names=FALSE)
tx2gene <- read.csv(file.path(BASE, "stringtie_rf/tx2gene_real_genes.csv"))
pheno <- read.csv(file.path(BASE, "ballgown_v3/phenotype.csv"))
pheno$tissue_type <- ifelse(grepl("N$", pheno$sample_id), "Normal", "Tumor")

cat(sprintf("  Transcripts: %d, Samples: %d\n", nrow(tx_counts), ncol(tx_counts)))
cat(sprintf("  Known gene symbols: %d unique\n", length(unique(tx2gene$gene_name[!grepl("^MSTRG", tx2gene$gene_name)]))))

# ── Filter: known transcripts only (XM_, NM_, XR_) ────────────────
cat("\n[2/5] Filtering to known RefSeq transcripts...\n")
known_tx <- grep("^[XN]M_|^XR_", rownames(tx_counts), value=TRUE)
cat(sprintf("  Known transcripts: %d\n", length(known_tx)))
tx_counts_known <- tx_counts[known_tx, ]

# Map to gene symbols
tx2gene_map <- setNames(tx2gene$gene_name, tx2gene$transcript_id)
gene_symbols <- tx2gene_map[rownames(tx_counts_known)]

# Filter: keep transcripts with ≥10 counts in ≥9 samples (25% of samples)
keep_tx <- rowSums(tx_counts_known >= 10) >= 9
tx_counts_filt <- tx_counts_known[keep_tx, ]
gene_symbols_filt <- gene_symbols[keep_tx]
cat(sprintf("  Expressed known transcripts: %d (%d unique genes)\n", 
            nrow(tx_counts_filt), length(unique(gene_symbols_filt))))

# ── Isoform-level PCA ──────────────────────────────────────────────
cat("\n[3/5] Isoform-level PCA (top 500 variable transcripts)...\n")

# DESeq2 VST on all samples (blind)
pheno$timepoint <- factor(pheno$timepoint, levels=c("Day0","Day2","Day6","Recovered"))
pheno$dog_id <- factor(pheno$dog_id)
rownames(pheno) <- pheno$sample_id

dds_tx <- DESeqDataSetFromMatrix(
  countData = tx_counts_filt[, pheno$sample_id],
  colData = pheno,
  design = ~ 1
)

# Remove genes with all-zero rows (shouldn't exist after filter, but safe)
dds_tx <- dds_tx[rowSums(counts(dds_tx)) > 0, ]
dds_tx <- estimateSizeFactors(dds_tx)

cat(sprintf("  VST matrix: %d transcripts x %d samples\n", nrow(dds_tx), ncol(dds_tx)))

vsd_tx <- vst(dds_tx, blind=TRUE)
log_tx <- assay(vsd_tx)

# Top 500 most variable transcripts
rv_tx <- rowVars(log_tx)
select_500 <- head(order(rv_tx, decreasing=TRUE), 500)

pca_tx <- prcomp(t(log_tx[select_500, ]), center=TRUE, scale.=FALSE)
var_exp <- pca_tx$sdev^2 / sum(pca_tx$sdev^2) * 100

cat(sprintf("  Isoform PCA: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%\n", 
            var_exp[1], var_exp[2], var_exp[3]))

pca_df <- as.data.frame(pca_tx$x)
pca_df$sample <- rownames(pca_df)
pca_df$dog <- pheno$dog_id
pca_df$timepoint <- pheno$timepoint
pca_df$tissue <- pheno$tissue_type

# ── PLOT: Isoform PCA ──────────────────────────────────────────────
cat("\n[4/5] Plotting isoform PCA...\n")

tp_colors <- c(Day0="#ef4444", Day2="#f97316", Day6="#eab308", Recovered="#22c55e")
dog_colors <- c(T1="#6366f1", T2="#8b5cf6", T3="#a855f7", T4="#d946ef",
                T5="#ec4899", T6="#f43f5e", T7="#14b8a6", T8="#0ea5e9", T9="#84cc16")

p1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=timepoint, shape=tissue)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=dog), size=2.8, max.overlaps=30, show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  scale_shape_manual(values=c(Tumor=16, Normal=17)) +
  stat_ellipse(aes(group=tissue), type="norm", level=0.95, linetype="dashed", alpha=0.3) +
  labs(title="CTVT: Isoform-Level PCA (Known RefSeq Transcripts)",
       subtitle=sprintf("%d known transcripts, 500 most variable • RefSeq: XM_/NM_/XR_", nrow(tx_counts_filt)),
       x=sprintf("PC1 (%.1f%%)", var_exp[1]),
       y=sprintf("PC2 (%.1f%%)", var_exp[2]),
       color="Timepoint", shape="Tissue") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text(hjust=0.5, color="grey40"))

ggsave(file.path(OUTDIR, "isoform_pca_timepoint.png"), p1, width=11, height=8, dpi=150)

p2 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=dog, shape=tissue)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=timepoint), size=2.5, max.overlaps=30, show.legend=FALSE) +
  scale_color_manual(values=dog_colors) +
  scale_shape_manual(values=c(Tumor=16, Normal=17)) +
  labs(title="CTVT: Isoform PCA by Dog",
       x=sprintf("PC1 (%.1f%%)", var_exp[1]),
       y=sprintf("PC2 (%.1f%%)", var_exp[2]),
       color="Dog", shape="Tissue") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5))

ggsave(file.path(OUTDIR, "isoform_pca_by_dog.png"), p2, width=11, height=8, dpi=150)

p3 <- ggplot(pca_df, aes(x=PC3, y=PC4, color=timepoint, shape=tissue)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(aes(label=paste0(dog, ifelse(tissue=="Normal","(N)",""))), 
                  size=2.5, max.overlaps=30, show.legend=FALSE) +
  scale_color_manual(values=tp_colors) +
  scale_shape_manual(values=c(Tumor=16, Normal=17)) +
  labs(title="CTVT: Isoform PCA — PC3 vs PC4",
       x=sprintf("PC3 (%.1f%%)", var_exp[3]),
       y=sprintf("PC4 (%.1f%%)", var_exp[4])) +
  theme_minimal(base_size=13)

ggsave(file.path(OUTDIR, "isoform_pca_pc3pc4.png"), p3, width=11, height=8, dpi=150)

# ── Gene-level proportion DTU (chi-squared test) ──────────────────
cat("\n[5/5] Running gene-level proportion DTU test...\n")
cat("  (chi-squared on per-gene isoform proportions, Tumor vs Normal)\n")

# For each gene with ≥2 expressed isoforms:
# Compare isoform proportion distributions between Tumor and Normal
gene_list <- unique(gene_symbols_filt[rownames(log_tx)])

dtu_results <- data.frame()
g1_samples <- which(pheno$tissue_type == "Tumor")
g2_samples <- which(pheno$tissue_type == "Normal")

for (g in gene_list) {
  tx_for_gene <- rownames(log_tx)[gene_symbols_filt[rownames(log_tx)] == g]
  if (length(tx_for_gene) < 2) next
  
  # Sum counts per transcript per group
  g1_counts <- rowSums(counts(dds_tx)[tx_for_gene, g1_samples, drop=FALSE])
  g2_counts <- rowSums(counts(dds_tx)[tx_for_gene, g2_samples, drop=FALSE])
  
  total <- sum(g1_counts) + sum(g2_counts)
  if (total < 20) next  # too few counts
  
  # Chi-squared test on the 2×n transcript contingency table
  tbl <- cbind(g1_counts, g2_counts)
  # Remove rows with zeros in both groups
  tbl <- tbl[rowSums(tbl) > 0, , drop=FALSE]
  if (nrow(tbl) < 2) next
  
  # Chi-squared test
  csq <- suppressWarnings(chisq.test(tbl))
  
  dtu_results <- rbind(dtu_results, data.frame(
    gene = g,
    n_transcripts = nrow(tbl),
    total_counts = total,
    chi_sq = csq$statistic,
    df = csq$parameter,
    pvalue = csq$p.value,
    stringsAsFactors = FALSE
  ))
}

# BH correction
dtu_results$padj <- p.adjust(dtu_results$pvalue, method="BH")
dtu_results <- dtu_results[order(dtu_results$pvalue), ]

sig_dtu <- sum(dtu_results$padj < 0.05, na.rm=TRUE)
cat(sprintf("  Genes tested: %d, Significant DTU (padj<0.05): %d\n", 
            nrow(dtu_results), sig_dtu))

write.csv(dtu_results, file.path(OUTDIR, "isoform_dtu_genes.csv"), row.names=FALSE)

# Top DTU genes
cat("\n  Top 15 DTU genes:\n")
for (i in 1:min(15, nrow(dtu_results))) {
  r <- dtu_results[i, ]
  cat(sprintf("  %2d. %-25s χ²=%.1f  p=%.2e  padj=%.2e\n", 
              i, r$gene, r$chi_sq, r$pvalue, r$padj))
}

# ── Heatmap: Top DTU gene isoforms ──────────────────────────────────
cat("\n  Generating DTU heatmap...\n")
top_dtu_genes <- head(dtu_results$gene[dtu_results$padj < 0.05], 50)
if (length(top_dtu_genes) > 0) {
  top_dtu_tx <- rownames(log_tx)[gene_symbols_filt[rownames(log_tx)] %in% top_dtu_genes]
  if (length(top_dtu_tx) > 0 && length(top_dtu_tx) <= 500) {
    hm_data <- log_tx[top_dtu_tx, ]
    hm_scaled <- t(scale(t(hm_data)))
    hm_scaled[hm_scaled > 3] <- 3
    hm_scaled[hm_scaled < -3] <- -3
    
    annotation_col <- data.frame(
      Timepoint = pheno$timepoint,
      Tissue = pheno$tissue_type,
      row.names = colnames(hm_scaled)
    )
    ann_colors <- list(
      Timepoint = tp_colors,
      Tissue = c(Tumor="#ef4444", Normal="#22c55e")
    )
    
    # Labels: transcript_id (gene_symbol)
    tx_labels <- paste0(substr(top_dtu_tx, 1, 18), " (", 
                        gene_symbols_filt[top_dtu_tx], ")")
    
    png(file.path(OUTDIR, "isoform_dtu_heatmap.png"),
        width=14, height=max(12, length(top_dtu_tx)*0.2), units="in", res=150)
    pheatmap(hm_scaled,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             labels_row = tx_labels,
             show_rownames = nrow(hm_scaled) <= 100,
             show_colnames = TRUE,
             fontsize_row = 5,
             fontsize_col = 7,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             main = "CTVT: Top DTU Isoforms (Tumor vs Normal)",
             border_color = NA)
    dev.off()
    cat(sprintf("  Heatmap: %d transcripts from top %d DTU genes\n", 
                length(top_dtu_tx), length(top_dtu_genes)))
  }
}

# ── Save everything ────────────────────────────────────────────────
write.csv(pca_df, file.path(OUTDIR, "isoform_pca_scores.csv"), row.names=FALSE)
write.csv(as.data.frame(log_tx), file.path(OUTDIR, "isoform_vst_expression.csv"))

cat("\n═══════════════════════════════════════════\n")
cat("  ✅ Isoform analysis complete\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat(sprintf("  Files: %s\n", paste(list.files(OUTDIR), collapse=", ")))

# ── Telegram Summary ────────────────────────────────────────────────
cat("\n── TELEGRAM SUMMARY ──\n")
cat("📊 CTVT Isoform-Level Analysis\n")
cat(sprintf("  Transcripts: %d RefSeq (filtered from %d known)\n", 
            nrow(tx_counts_filt), length(known_tx)))
cat(sprintf("  Unique genes: %d\n", length(gene_list)))
cat(sprintf("  PCA: PC1=%.1f%% (Tumor vs Normal), PC2=%.1f%%\n", var_exp[1], var_exp[2]))
cat(sprintf("  DTU (chi-sq): %d genes tested, %d significant\n", 
            nrow(dtu_results), sig_dtu))
cat(sprintf("  DTU p<0.05 uncorrected: %d genes\n",
            sum(dtu_results$pvalue < 0.05, na.rm=TRUE)))
