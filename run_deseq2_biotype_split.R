#!/usr/bin/env Rscript
# ====================================================================
# CTVT RNA-seq: DESeq2 Split by Biotype — Protein-coding & Known lncRNA
# Uses Salmon + tximport (no featureCounts needed)
# Design: ~ dog_id + timepoint (paired longitudinal, 9 dogs × 4 timepoints)
# ====================================================================

Sys.setenv(R_MAX_VSIZE = "250e9")
options(future.globals.maxSize = 10 * 1024^3)

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"
OUTLIERS <- c("T1B_S2", "T8C_S10", "T9N_S26")  # Old featureCounts-based flags (for comparison only)

dir.create(file.path(BASE, "deseq2_biotype_split"), showWarnings = FALSE, recursive = TRUE)

# ── 1. Sample metadata ─────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════\n")
cat("  CTVT DESeq2 — Biotype-Split Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

all_samples <- basename(list.dirs(QUANT_DIR, recursive = FALSE, full.names = FALSE))
keep_samples <- setdiff(all_samples, OUTLIERS)

sample_meta <- data.frame(sampleID = keep_samples, stringsAsFactors = FALSE)
sample_meta$timepoint <- factor(sapply(sample_meta$sampleID, function(s) {
  c(A="Day0",B="Day2",C="Day6",N="Recovered")[substr(s,3,3)]
}), levels = c("Day0", "Day2", "Day6", "Recovered"))

sample_meta$dog_id <- factor(sapply(sample_meta$sampleID, function(s) {
  substr(strsplit(s, "_")[[1]][1], 1, 2)
}))

cat(sprintf("  %d samples: Day0=%d, Day2=%d, Day6=%d, Recovered=%d\n",
  nrow(sample_meta),
  sum(sample_meta$timepoint == "Day0"),
  sum(sample_meta$timepoint == "Day2"),
  sum(sample_meta$timepoint == "Day6"),
  sum(sample_meta$timepoint == "Recovered")))

# ── 2. Define biotype runs ─────────────────────────────────────────
biotype_runs <- list(
  protein_coding = list(
    tx2gene_file = file.path(BASE, "tx2gene_protein_coding.csv"),
    label = "Protein-Coding",
    outfile = "protein_coding"
  ),
  lncrna_known = list(
    tx2gene_file = file.path(BASE, "tx2gene_lncrna_known.csv"),
    label = "Known lncRNA",
    outfile = "lncrna_known"
  )
)

contrasts <- list(
  Day2_vs_Day0 = c("Day2", "Day0"),
  Day6_vs_Day0 = c("Day6", "Day0"),
  Day6_vs_Day2 = c("Day6", "Day2"),
  Recovered_vs_Day0 = c("Recovered", "Day0"),
  Recovered_vs_Day6 = c("Recovered", "Day6")
)

# ── 3. Run DESeq2 for each biotype ─────────────────────────────────
for (bt_name in names(biotype_runs)) {
  cfg <- biotype_runs[[bt_name]]
  outdir <- file.path(BASE, "deseq2_biotype_split", cfg$outfile)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  cat(sprintf("\n\n%s %s %s\n", strrep("═", 20), cfg$label, strrep("═", 20)))
  
  # Load tx2gene
  tx2gene <- read_csv(cfg$tx2gene_file, show_col_types = FALSE)
  cat(sprintf("  tx2gene: %d transcripts → %d genes\n", 
    nrow(tx2gene), length(unique(tx2gene$GENEID))))
  
  # Import Salmon quant
  files <- file.path(QUANT_DIR, sample_meta$sampleID, "quant.sf")
  names(files) <- sample_meta$sampleID
  
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                   countsFromAbundance = "lengthScaledTPM")
  
  cat(sprintf("  Imported: %d genes x %d samples\n", nrow(txi$counts), ncol(txi$counts)))
  
  # Filter low counts
  keep <- rowSums(txi$counts >= 10) >= 4
  txi_filt <- lapply(txi, function(x) {
    if (is.matrix(x)) x[keep, , drop = FALSE] else x
  })
  cat(sprintf("  After filtering: %d genes\n", nrow(txi_filt$counts)))
  
  if (nrow(txi_filt$counts) < 100) {
    cat("  ⚠️ Too few genes, skipping\n")
    next
  }
  
  # DESeq2
  dds <- DESeqDataSetFromTximport(txi_filt, sample_meta, ~ dog_id + timepoint)
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract all contrasts
  all_results <- list()
  deg_summary <- data.frame()
  
  for (cname in names(contrasts)) {
    c1 <- contrasts[[cname]][1]
    c2 <- contrasts[[cname]][2]
    
    res <- results(dds, contrast = c("timepoint", c1, c2), alpha = 0.05)
    res_shrunk <- lfcShrink(dds, contrast = c("timepoint", c1, c2), 
                             res = res, type = "normal")
    res_df <- as.data.frame(res)
    res_shrunk_df <- as.data.frame(res_shrunk)
    res_df$shrunkenLog2FC <- res_shrunk_df$log2FoldChange
    res_df$shrunkenlfcSE <- res_shrunk_df$lfcSE
    res_df$gene <- rownames(res_df)
    
    n_deg <- sum(res$padj < 0.05, na.rm = TRUE)
    n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
    n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
    
    cat(sprintf("  %s: %d DEGs (%d up, %d down)\n", cname, n_deg, n_up, n_down))
    
    all_results[[cname]] <- res_df
    
    deg_summary <- rbind(deg_summary, data.frame(
      Contrast = cname, DE_genes = n_deg, Up = n_up, Down = n_down
    ))
    
    # Save significant DEGs with shrunken LFC
    sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
    sig <- sig[order(sig$padj), ]
    write.csv(sig, file.path(outdir, paste0(cname, "_DEGs.csv")), row.names = TRUE)
    
    # Save ALL results (for GSEA/ranking)
    all_res <- res_df[order(res_df$padj, na.last = TRUE), ]
    write.csv(all_res, file.path(outdir, paste0(cname, "_all_results.csv")), row.names = TRUE)
  }
  
  # Save summary
  write.csv(deg_summary, file.path(outdir, "DE_summary.csv"), row.names = FALSE)
  
  cat(sprintf("\n  ✅ %s results saved to %s\n", cfg$label, outdir))
}

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  All biotype-split DESeq2 runs complete!\n")
cat("═══════════════════════════════════════════════════════════════\n")