#!/usr/bin/env Rscript
# ====================================================================
# CTVT DRIMSeq v5 — Fixed filtering for isoform-level DTU
# Problem with v4: G-test underflow, everything significant
# Problem with DRIMSeq v3: filters too strict, killing real signals
# Fix: Use DEXSeq-like filtering, test ALL 5 tumor-timepoint contrasts
# ====================================================================

suppressPackageStartupMessages({
  library(DRIMSeq)
  library(stageR)
  library(ggplot2)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
OUTDIR <- file.path(BASE, "isoform_drimseq_v5")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════════\n")
cat("  CTVT DRIMSeq v5 — Isoform DTU Analysis\n")
cat("═══════════════════════════════════════════\n")

# ── Load data ──────────────────────────────────────────────────────
cat("\n[1/5] Loading data...\n")
counts_raw <- read.csv(file.path(BASE, "ballgown_v3/transcript_count_matrix_prepDE.csv"),
                       row.names=1, check.names=FALSE)
tx2gene <- read.csv(file.path(BASE, "stringtie_rf/tx2gene.csv"))
pheno <- read.csv(file.path(BASE, "ballgown_v3/phenotype.csv"))

samps <- data.frame(sample_id=pheno$sample_id, group=pheno$timepoint, 
                    dog=pheno$dog_id, stringsAsFactors=FALSE)
samps$group <- factor(samps$group, levels=c("Day0","Day2","Day6","Recovered"))
samps$dog <- factor(samps$dog)

counts_mat <- as.matrix(counts_raw[, as.character(samps$sample_id)])
tx2gene_sub <- tx2gene[match(rownames(counts_mat), tx2gene$transcript_id), ]

counts_df <- data.frame(gene_id=tx2gene_sub$gene_id,
                        feature_id=tx2gene_sub$transcript_id,
                        counts_mat, check.names=FALSE)
counts_df <- counts_df[!is.na(counts_df$gene_id), ]

cat(sprintf("  Loaded: %d genes × %d transcripts × %d samples\n",
            length(unique(counts_df$gene_id)), nrow(counts_df), nrow(samps)))

# ── GLOBAL FILTER: relaxed but rational ────────────────────────────
cat("\n[2/5] Global filtering...\n")

# For each gene, keep if it has ≥2 expressed isoforms
# "Expressed" = ≥5 counts in ≥12 samples (1/3 of all samples)
d <- dmDSdata(counts=counts_df, samples=samps)
d <- dmFilter(d, 
  min_samps_feature_expr = 18,   # expressed in at least half the samples
  min_feature_expr = 10,          # at least 10 counts
  min_feature_prop = 0.01,       # at least 1% of gene's reads
  min_samps_feature_prop = 18,
  min_samps_gene_expr = 18,      # gene expressed in at least 18 samples
  min_gene_expr = 10
)

cts <- counts(d); sams_global <- DRIMSeq::samples(d)
cat(sprintf("  After global filter: %d genes (%d transcripts)\n", 
            length(unique(cts$gene_id)), nrow(cts)))

# ── PER-CONTRAST ANALYSIS ──────────────────────────────────────────
cat("\n[3/5] Running per-contrast DRIMSeq...\n")

contrasts <- list(
  c("Day0", "Day2", "Day0_vs_Day2"),
  c("Day0", "Day6", "Day0_vs_Day6"),
  c("Day0", "Recovered", "Day0_vs_Recovered"),
  c("Day2", "Day6", "Day2_vs_Day6"),
  c("Day6", "Recovered", "Day6_vs_Recovered")
)

results_summary <- data.frame()
all_stager <- list()

for (i in seq_along(contrasts)) {
  g1 <- contrasts[[i]][1]
  g2 <- contrasts[[i]][2]
  cn <- contrasts[[i]][3]
  
  cat(sprintf("\n── %s (%s vs %s) ──\n", cn, g1, g2))
  
  # Subset samples for this contrast
  keep <- sams_global$group %in% c(g1, g2)
  sams_sub <- sams_global[keep, , drop=FALSE]
  sams_sub$group <- droplevels(sams_sub$group)
  
  n1 <- sum(sams_sub$group == g1)
  n2 <- sum(sams_sub$group == g2)
  cat(sprintf("  Samples: %s=%d, %s=%d\n", g1, n1, g2, n2))
  
  cts_sub <- cts[, c("gene_id", "feature_id", as.character(sams_sub$sample_id))]
  
  # Per-contrast filter: keep genes with ≥2 isoforms expressed
  # in at least half the samples of each group
  min_samps_contrast <- max(3, floor(min(n1, n2) / 2))
  
  d_sub <- dmDSdata(counts=cts_sub, samples=sams_sub)
  d_sub <- dmFilter(d_sub, 
    min_samps_feature_expr = min_samps_contrast,
    min_feature_expr = 5,
    min_feature_prop = 0.01,
    min_samps_feature_prop = min_samps_contrast,
    min_samps_gene_expr = min_samps_contrast * 2,  # total across both groups
    min_gene_expr = 10
  )
  
  ng <- length(unique(counts(d_sub)$gene_id))
  nt <- nrow(counts(d_sub))
  cat(sprintf("  After contrast filter: %d genes (%d transcripts)\n", ng, nt))
  
  if (ng < 10) {
    cat("  ⚠ Too few genes, skipping\n")
    next
  }
  
  # Design matrix
  design_sub <- model.matrix(~ group, data=DRIMSeq::samples(d_sub))
  coef_name <- setdiff(colnames(design_sub), "(Intercept)")
  cat(sprintf("  Coef: %s\n", coef_name))
  
  # Precision estimation (with retry)
  set.seed(42)
  tryCatch({
    d_sub <- dmPrecision(d_sub, design=design_sub, 
                         verbose=0, add_uniform=TRUE)
    d_sub <- dmFit(d_sub, design=design_sub, verbose=0)
    d_sub <- dmTest(d_sub, coef=coef_name, verbose=0)
  }, error=function(e) {
    cat(sprintf("  ⚠ DRIMSeq error: %s\n", e$message))
    d_sub <<- NULL
  })
  
  if (is.null(d_sub)) next
  
  # Extract results
  res_g <- results(d_sub)
  res_t <- results(d_sub, level="feature")
  
  sig_g_drim <- sum(res_g$adj_pvalue < 0.05, na.rm=TRUE)
  sig_t_drim <- sum(res_t$adj_pvalue < 0.05, na.rm=TRUE)
  cat(sprintf("  DRIMSeq raw: %d genes, %d tx (adj_p < 0.05)\n", sig_g_drim, sig_t_drim))
  
  # Save DRIMSeq results
  write.csv(res_g, file.path(OUTDIR, paste0(cn, "_genes.csv")), row.names=FALSE)
  write.csv(res_t, file.path(OUTDIR, paste0(cn, "_transcripts.csv")), row.names=FALSE)
  
  # ── stageR ──
  na_ok <- function(x) { x[is.na(x)] <- 1; x }
  pScreen <- na_ok(res_g$pvalue)
  names(pScreen) <- res_g$gene_id
  
  pConf <- matrix(na_ok(res_t$pvalue), ncol=1)
  rownames(pConf) <- res_t$feature_id
  colnames(pConf) <- "tx"
  
  tx2g <- unique(data.frame(
    feature_id=res_t$feature_id,
    gene_id=res_t$gene_id,
    stringsAsFactors=FALSE
  ))
  
  sr <- tryCatch({
    stageRTx(pScreen=pScreen, pConfirmation=pConf,
             pScreenAdjusted=FALSE, tx2gene=tx2g)
  }, error=function(e) {
    cat(sprintf("  ⚠ stageR init error: %s\n", e$message))
    NULL
  })
  
  if (is.null(sr)) next
  
  sr <- stageWiseAdjustment(sr, method="dtu", alpha=0.05, allowNA=TRUE)
  sp <- getAdjustedPValues(sr, order=FALSE, onlySignificantGenes=FALSE)
  
  sig_sr_g <- sum(sp$gene < 0.05, na.rm=TRUE)
  sig_sr_t <- sum(sp$transcript < 0.05, na.rm=TRUE)
  cat(sprintf("  stageR (gene-level corrected): %d genes, %d tx\n", sig_sr_g, sig_sr_t))
  
  write.csv(sp, file.path(OUTDIR, paste0(cn, "_stager.csv")), row.names=FALSE)
  all_stager[[cn]] <- sp
  
  results_summary <- rbind(results_summary, data.frame(
    Contrast = cn,
    Genes_Tested = ng,
    DRIMSeq_Genes = sig_g_drim,
    DRIMSeq_Tx = sig_t_drim,
    stageR_Genes = sig_sr_g,
    stageR_Tx = sig_sr_t,
    stringsAsFactors = FALSE
  ))
}

# ── Print summary ──────────────────────────────────────────────────
cat("\n[4/5] ═══════════════════════════════════════════\n")
cat("  RESULTS SUMMARY\n")
cat("═══════════════════════════════════════════\n")
print(results_summary, row.names=FALSE)
write.csv(results_summary, file.path(OUTDIR, "summary.csv"), row.names=FALSE)

# ── Per-contrast top genes ────────────────────────────────────────
cat("\n[5/5] Top DTU genes per contrast:\n")
for (cn in names(all_stager)) {
  sp <- all_stager[[cn]]
  sp_sig <- sp[sp$gene < 0.05 & !is.na(sp$gene), ]
  top5 <- head(sp_sig[order(sp_sig$gene), ], 5)
  cat(sprintf("\n%s: %d significant genes\n", cn, nrow(sp_sig)))
  if (nrow(top5) > 0) {
    for (j in 1:nrow(top5)) {
      cat(sprintf("  %s: gene_p=%.2e, tx_p=%.2e\n",
                  top5$geneID[j], top5$gene[j], top5$transcript[j]))
    }
  }
}

# ── Venn-style overlap ────────────────────────────────────────────
cat("\n── DTU gene overlap ──\n")
dtu_genes <- list()
for (cn in names(all_stager)) {
  sp <- all_stager[[cn]]
  dtu_genes[[cn]] <- unique(sp$geneID[sp$gene < 0.05 & !is.na(sp$gene)])
  cat(sprintf("  %s: %d genes\n", cn, length(dtu_genes[[cn]])))
}

# Pairwise
for (i in 1:(length(dtu_genes)-1)) {
  for (j in (i+1):length(dtu_genes)) {
    o <- intersect(dtu_genes[[i]], dtu_genes[[j]])
    if (length(o) > 0) {
      cat(sprintf("  %s ∩ %s: %d\n", names(dtu_genes)[i], names(dtu_genes)[j], length(o)))
    }
  }
}

# ── Save DTU gene lists ──
for (cn in names(dtu_genes)) {
  writeLines(dtu_genes[[cn]], file.path(OUTDIR, paste0(cn, "_dtu_genes.txt")))
}

cat(sprintf("\n✅ Results saved to: %s/\n", OUTDIR))
cat(sprintf("   Files: %s\n", paste(list.files(OUTDIR), collapse=", ")))
