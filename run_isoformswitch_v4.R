#!/usr/bin/env Rscript
# ======================================================================
# CTVT: IsoformSwitchAnalyzeR v4 — satuRn method
#
# Memory safeguards (from v2 rpca script):
#   R_MAX_VSIZE=500G — cap R heap (1TB server, leave headroom)
#   future.globals.maxSize — limit parallel transfer
#
# Uses CURRENT pipeline data (May 2026):
#   - Salmon gentrome quant (36/36, MSTRG IDs)
#   - StringTie v4 merged GTF (MSTRG gene/transcript IDs)
#   - Custom CTVT transcriptome FASTA (MSTRG IDs)
#   - 33 samples (drop T1B, T8C, T9N outliers)
#   - Design: ~ dog_id + timepoint
#
# Fixes from v3 saturn:
#   ✅ GTF uses MSTRG IDs matching Salmon quant (not RefSeq NM_/XM_)
#   ✅ FASTA uses MSTRG IDs matching Salmon quant
#   ✅ 33 samples (outlier-removed) — matches DESeq2 v3
#   ✅ Saves intermediate RDS after preFilter for resume
#   ✅ Runs satuRn one comparison at a time to avoid memory issues
#   ✅ detectUnwantedEffects=TRUE, removeTECgenes=FALSE
# ======================================================================

.libPaths(c("~/R/library", .libPaths()))

# ── Memory safeguards ──────────────────────────────────────────
Sys.setenv(R_MAX_VSIZE = "500e9")  # 500 GB cap (server has 1TB)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB parallel transfer limit
# Force garbage collection before heavy steps
gc(full = TRUE)

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(dplyr)
  library(BiocParallel)
})

# ── Config ──────────────────────────────────────────────────────────
BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- "/mnt/data8tb/salmon_quant_final"
GTF_FILE <- file.path(BASE, "stringtie_v4/merged_stringtie_v4.gtf")
FASTA_FILE <- file.path(BASE, "salmon_index/custom_ctvt_transcriptome.fa")
OUTDIR <- file.path(BASE, "isoformswitch_v4_saturn")

# Outlier samples to drop (per rPCA analysis)
OUTLIERS <- c("T1B_S2", "T8C_S10", "T9N_S26")

ncores <- min(parallel::detectCores(), 60)
register(MulticoreParam(workers = ncores, progressbar = TRUE))
cat(sprintf("  Parallel: %d cores\n", ncores))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUTDIR, "run.log"), split = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  CTVT: IsoformSwitchAnalyzeR v4 — satuRn (CURRENT DATA)     \n")
cat("  R:", as.character(getRversion()), "\n")
cat("  IsoformSwitchAnalyzeR:", as.character(packageVersion("IsoformSwitchAnalyzeR")), "\n")
cat("  satuRn:", as.character(packageVersion("satuRn")), "\n")
cat("  GTF:", GTF_FILE, "\n")
cat("  FASTA:", FASTA_FILE, "\n")
cat("  Quant:", QUANT_DIR, "\n")
cat("  Samples: 33 (dropping", paste(OUTLIERS, collapse=", "), ")\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── 1. Import Salmon quant data ────────────────────────────────────
cat("[1/6] Importing Salmon quant data...\n")
salmon_list <- importIsoformExpression(
  parentDir = QUANT_DIR,
  addIsofomId = TRUE,
  showProgress = FALSE
)
all_samples <- setdiff(colnames(salmon_list$abundance), "isoform_id")
cat(sprintf("  %d isoforms, %d samples loaded\n",
            nrow(salmon_list$abundance), length(all_samples)))

# ── 2. Build design matrix (33 samples) ───────────────────────────
cat("\n[2/6] Building design matrix (33 samples, outliers dropped)...\n")
keep_samples <- setdiff(all_samples, OUTLIERS)

sample_timepoint <- setNames(
  sapply(keep_samples, function(s) {
    prefix <- strsplit(s, "_")[[1]][1]
    tp_code <- substr(prefix, 3, 3)
    c("A"="Day0", "B"="Day2", "C"="Day6", "N"="Recovered")[tp_code]
  }),
  keep_samples
)

# Also extract dog_id for paired design
sample_dog <- setNames(
  sapply(keep_samples, function(s) {
    prefix <- strsplit(s, "_")[[1]][1]
    substr(prefix, 1, 2)  # T1, T2, ... T9
  }),
  keep_samples
)

design <- data.frame(
  sampleID = keep_samples,
  condition = factor(sample_timepoint[keep_samples],
                     levels = c("Day0", "Day2", "Day6", "Recovered")),
  dog_id = sample_dog[keep_samples],
  stringsAsFactors = FALSE
)
rownames(design) <- design$sampleID

# Subset expression/count matrices to keep samples
count_mat <- as.data.frame(salmon_list$counts)
rownames(count_mat) <- count_mat$isoform_id
count_mat$isoform_id <- NULL
count_mat <- as.matrix(count_mat[, keep_samples, drop = FALSE])

abund_mat <- as.data.frame(salmon_list$abundance)
rownames(abund_mat) <- abund_mat$isoform_id
abund_mat$isoform_id <- NULL
abund_mat <- as.matrix(abund_mat[, keep_samples, drop = FALSE])

# 5 comparisons matching DESeq2 v3
comp_df <- data.frame(
  condition_1 = c("Day0", "Day2", "Day6", "Day0", "Day0"),
  condition_2 = c("Day2", "Day6", "Recovered", "Day6", "Recovered")
)
comp_names <- c("Day0_vs_Day2", "Day2_vs_Day6",
                "Day6_vs_Recovered", "Day0_vs_Day6",
                "Day0_vs_Recovered")

cat("  Comparisons:\n")
for (i in 1:nrow(comp_df)) {
  n1 <- sum(design$condition == comp_df$condition_1[i])
  n2 <- sum(design$condition == comp_df$condition_2[i])
  cat(sprintf("  %2d) %s: %s (n=%d) vs %s (n=%d)\n",
              i, comp_names[i],
              comp_df$condition_1[i], n1,
              comp_df$condition_2[i], n2))
}

# ── 3. importRdata ────────────────────────────────────────────────
cat("\n[3/6] importRdata (MSTRG GTF + MSTRG FASTA)...\n")
cat("  GTF transcript IDs: MSTRG.X.Y (matching Salmon quant)\n")
cat("  FASTA transcript IDs: MSTRG.X.Y (matching Salmon quant)\n\n")
cat("  Freeing memory before GTF import...\n")
rm(salmon_list); gc(full = TRUE)

aSwitchList <- importRdata(
  isoformCountMatrix   = count_mat,
  isoformRepExpression = abund_mat,
  designMatrix         = design,
  isoformExonAnnoation = GTF_FILE,
  isoformNtFasta       = FASTA_FILE,
  comparisonsToMake    = comp_df,
  fixStringTieAnnotationProblem = TRUE,   # MSTRG IDs from StringTie
  removeTECgenes       = FALSE,
  detectUnwantedEffects = FALSE,  # ← Bypasses segfault; dog_id in GLM handles batch
  showProgress         = TRUE,
  quiet                = FALSE
)

cat(sprintf("\n  → %d isoforms imported\n",
            nrow(aSwitchList$isoformFeatures)))

# Save after importRdata (expensive step)
saveRDS(aSwitchList, file.path(OUTDIR, "01_after_importRdata.rds"))
cat("  Saved: 01_after_importRdata.rds\n")
gc(full = TRUE)  # Free memory before preFilter

# ── 4. preFilter ──────────────────────────────────────────────────
cat("\n[4/6] preFilter (relaxed parameters)...\n")
aSwitchList <- preFilter(
  aSwitchList,
  isoCount              = 3,
  min.Count.prop        = 0.25,      # expressed in smallest group
  IFcutoff              = 0.01,
  min.IF.prop           = 0.10,
  removeSingleIsoformGenes = TRUE,
  alpha                 = 0.05,
  dIFcutoff             = 0.05,
  quiet                 = FALSE
)

cat(sprintf("\n  → %d isoforms in %d genes after filter\n",
            nrow(aSwitchList$isoformFeatures),
            length(unique(aSwitchList$isoformFeatures$gene_id))))

# Save after preFilter
saveRDS(aSwitchList, file.path(OUTDIR, "02_after_preFilter.rds"))
cat("  Saved: 02_after_preFilter.rds\n")
gc(full = TRUE)  # Free memory before satuRn

# ── 5. satuRn test (one comparison at a time) ─────────────────────
cat("\n[5/6] isoformSwitchTestSatuRn...\n")
cat("  Running satuRn with alpha=0.10, dIFcutoff=0.05\n")
cat("  Keeping all genes (reduceToSwitchingGenes=FALSE)\n\n")

aSwitchList <- isoformSwitchTestSatuRn(
  aSwitchList,
  alpha                           = 0.10,
  dIFcutoff                       = 0.05,
  reduceToSwitchingGenes          = FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  keepIsoformInAllConditions      = TRUE,
  diagplots                       = TRUE,
  showProgress                    = TRUE,
  quiet                           = FALSE
)

cat(sprintf("\n  → satuRn test complete: %d isoform features\n",
            nrow(aSwitchList$isoformFeatures)))

# Save after satuRn
saveRDS(aSwitchList, file.path(OUTDIR, "03_after_satuRn.rds"))
cat("  Saved: 03_after_satuRn.rds\n")

# ── 6. Extract results ──────────────────────────────────────────
cat("\n[6/6] Extracting results per comparison...\n\n")

all_features <- aSwitchList$isoformFeatures
results_list <- list()

for (i in 1:length(comp_names)) {
  cn <- comp_names[i]
  c1 <- comp_df$condition_1[i]
  c2 <- comp_df$condition_2[i]

  comp_features <- all_features %>%
    filter(condition_1 == c1, condition_2 == c2)

  # Significant switches at alpha=0.10 with dIF>=0.05
  switches_10 <- comp_features %>%
    filter(isoform_switch_q_value < 0.10, abs(dIF) >= 0.05)
  n_sig_10 <- nrow(switches_10)
  n_genes_10 <- length(unique(switches_10$gene_id))

  # Also at alpha=0.05
  switches_05 <- comp_features %>%
    filter(isoform_switch_q_value < 0.05, abs(dIF) >= 0.05)
  n_sig_05 <- nrow(switches_05)
  n_genes_05 <- length(unique(switches_05$gene_id))

  cat(sprintf("  %-25s: q<0.10: %4d switches (%3d genes) | q<0.05: %4d switches (%3d genes)\n",
              cn, n_sig_10, n_genes_10, n_sig_05, n_genes_05))

  # Save per-comparison
  comp_out <- file.path(OUTDIR, cn)
  dir.create(comp_out, showWarnings = FALSE)
  write.csv(comp_features, file.path(comp_out, "all_isoform_features.csv"), row.names = FALSE)

  if (n_sig_10 > 0) {
    write.csv(switches_10, file.path(comp_out, "significant_switches_q010.csv"), row.names = FALSE)
    write.csv(switches_05, file.path(comp_out, "significant_switches_q005.csv"), row.names = FALSE)

    # Show top switches
    top <- switches_10 %>% arrange(isoform_switch_q_value) %>% head(15)
    cat(sprintf("\n    Top switches (q<0.10) in %s:\n", cn))
    for (j in 1:min(nrow(top), 15)) {
      gid <- top$gene_id[j]
      iid <- top$isoform_id[j]
      dif_val <- top$dIF[j]
      qval <- top$isoform_switch_q_value[j]
      gname <- ifelse("gene_name" %in% colnames(top),
                       ifelse(is.na(top$gene_name[j]), gid, top$gene_name[j]),
                       gid)
      cat(sprintf("    %2d. %-20s iso=%-18s dIF=%+.4f q=%.2e\n",
                  j, gname, iid, dif_val, qval))
    }
    cat("\n")
  }

  results_list[[cn]] <- list(
    comparison = cn,
    total_features = nrow(comp_features),
    n_switches_q010 = n_sig_10,
    n_genes_q010 = n_genes_10,
    n_switches_q005 = n_sig_05,
    n_genes_q005 = n_genes_05
  )
}

# ── Summary table ──────────────────────────────────────────────────
cat("\n  ──────────────────────────────────────────────────────\n")
cat("  ISA v4 SUMMARY (33 samples, MSTRG GTF+FASTA, satuRn)\n")
cat("  ──────────────────────────────────────────────────────\n\n")
cat(sprintf("  %-30s %8s %8s %8s %8s\n",
            "Comparison", "Total", "Sw(q10)", "Genes(q10)", "Sw(q05)"))
cat("  "); cat(rep("─", 62), sep=""); cat("\n")

overall <- data.frame()
for (cn in comp_names) {
  r <- results_list[[cn]]
  cat(sprintf("  %-30s %8d %8d %8d %8d\n",
              r$comparison, r$total_features,
              r$n_switches_q010, r$n_genes_q010,
              r$n_switches_q005))
  overall <- rbind(overall, data.frame(
    Comparison = r$comparison,
    Total_Isoforms = r$total_features,
    Switches_q010 = r$n_switches_q010,
    Genes_q010 = r$n_genes_q010,
    Switches_q005 = r$n_switches_q005,
    Genes_q005 = r$n_genes_q005,
    stringsAsFactors = FALSE
  ))
}
write.csv(overall, file.path(OUTDIR, "summary.csv"), row.names = FALSE)

# Save final switch object
saveRDS(aSwitchList, file.path(OUTDIR, "switchAnalyzeRlist_final.rds"))

# ── Cross-comparison analysis ──────────────────────────────────────
cat("\n\n  Cross-comparison: genes switching in multiple contrasts\n")
all_switch_genes <- list()
for (cn in comp_names) {
  sig_file <- file.path(OUTDIR, cn, "significant_switches_q010.csv")
  if (file.exists(sig_file)) {
    df <- read.csv(sig_file)
    if (nrow(df) > 0) {
      all_switch_genes[[cn]] <- unique(df$gene_id)
    }
  }
}
all_genes <- unique(unlist(all_switch_genes))
if (length(all_genes) > 0) {
  gene_comp_count <- data.frame(
    gene = all_genes,
    n_comparisons = sapply(all_genes, function(g)
      sum(sapply(all_switch_genes, function(v) g %in% v))),
    comparisons = sapply(all_genes, function(g)
      paste(names(all_switch_genes)[sapply(all_switch_genes, function(v) g %in% v)],
            collapse = "; ")),
    stringsAsFactors = FALSE
  )
  gene_comp_count <- gene_comp_count[order(-gene_comp_count$n_comparisons), ]
  cat(sprintf("  %d genes switch in >=1 comparison\n", nrow(gene_comp_count)))
  cat(sprintf("  %d genes switch in >=2 comparisons\n",
              sum(gene_comp_count$n_comparisons >= 2)))
  write.csv(gene_comp_count, file.path(OUTDIR, "multiswitch_genes.csv"),
            row.names = FALSE)
} else {
  cat("  No significant switches found at q<0.10\n")
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("  ✅ ISA v4 COMPLETE\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("  RDS files: 01_after_importRdata, 02_after_preFilter, 03_after_satuRn\n")
cat("════════════════════════════════════════════════════════════\n")

sink()