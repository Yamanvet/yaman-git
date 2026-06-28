#!/usr/bin/env Rscript
# ======================================================================
# CTVT: IsoformSwitchAnalyzeR v2.12.0 — 5 Pairwise Comparisons
# Optimized: single GTF import with all 5 comparisons
# ======================================================================

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(dplyr)
  library(BiocParallel)
})

# Use 60 cores
ncores <- 60
register(MulticoreParam(workers = ncores, progressbar = TRUE))
cat(sprintf("  Parallel: %d cores\n", ncores))

BASE <- "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR <- file.path(BASE, "salmon_quant")
GTF_FILE <- "/mnt/data8tb/CTVT_Raw_Data_root_backup/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gtf"
FASTA_FILE <- file.path(BASE, "salmon_index/canine_transcriptome.fa")
OUTDIR <- file.path(BASE, "isoformswitch_results_v2")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("════════════════════════════════════════════════════════════╗\n")
cat("  CTVT: IsoformSwitchAnalyzeR v2.12.0 — 5 Comparisons    ║\n")
cat("  R:", as.character(getRversion()), "     BioC:", as.character(packageVersion("BiocGenerics")), "  ║\n")
cat("════════════════════════════════════════════════════════════╝\n\n")

# ── 1. Import all Salmon data (once) ──────────────────────────────
cat("[1/5] Importing all Salmon quant data...\n")
salmon_list <- importIsoformExpression(
  parentDir = QUANT_DIR,
  addIsofomId = TRUE,
  showProgress = FALSE
)
all_samples <- setdiff(colnames(salmon_list$abundance), "isoform_id")
cat(sprintf("  %d isoforms, %d samples loaded\n", nrow(salmon_list$abundance), length(all_samples)))

# Map timepoints
sample_timepoint <- setNames(
  sapply(all_samples, function(s) {
    prefix <- strsplit(s, "_")[[1]][1]
    tp_code <- substr(prefix, 3, 3)
    c("A"="Day0", "B"="Day2", "C"="Day6", "N"="Recovered")[tp_code]
  }),
  all_samples
)

# ── 2. Build full design matrix with ALL 36 samples ───────────────
cat("\n[2/5] Building design matrix (all 36 samples, 4 timepoints)...\n")
design <- data.frame(
  sampleID = all_samples,
  condition = sample_timepoint[all_samples],
  stringsAsFactors = FALSE
)
design$condition <- factor(design$condition, levels = c("Day0", "Day2", "Day6", "Recovered"))
rownames(design) <- design$sampleID

# Define all 5 pairwise comparisons
comp_df <- data.frame(
  condition_1 = c("Day0", "Day2", "Day6", "Day0", "Day0"),
  condition_2 = c("Day2", "Day6", "Recovered", "Day6", "Recovered")
)
comp_names <- c("Day0_vs_Day2", "Day2_vs_Day6", "Day6_vs_Recovered", "Day0_vs_Day6", "Day0_vs_Recovered")

cat("  Comparisons:\n")
for (i in 1:nrow(comp_df)) {
  n1 <- sum(design$condition == comp_df$condition_1[i])
  n2 <- sum(design$condition == comp_df$condition_2[i])
  cat(sprintf("  %2d) %s: %s (n=%d) vs %s (n=%d)\n",
              i, comp_names[i], comp_df$condition_1[i], n1, comp_df$condition_2[i], n2))
}

# ── 3. Single importRdata with all comparisons ─────────────────────
cat("\n[3/5] importRdata with RefSeq GTF (single pass)...\n")

count_mat <- as.data.frame(salmon_list$counts)
rownames(count_mat) <- count_mat$isoform_id
count_mat$isoform_id <- NULL
count_mat <- as.matrix(count_mat)

abund_mat <- as.data.frame(salmon_list$abundance)
rownames(abund_mat) <- abund_mat$isoform_id
abund_mat$isoform_id <- NULL
abund_mat <- as.matrix(abund_mat)

aSwitchList <- importRdata(
  isoformCountMatrix   = count_mat,
  isoformRepExpression = abund_mat,
  designMatrix         = design,
  isoformExonAnnoation = GTF_FILE,
  isoformNtFasta       = FASTA_FILE,
  comparisonsToMake    = comp_df,
  fixStringTieAnnotationProblem = FALSE,
  removeTECgenes = TRUE,
  detectUnwantedEffects = FALSE,   # Skip SVA — saves 10-20 min, known issue at scale
  showProgress = TRUE,
  quiet = FALSE
)

cat(sprintf("  → %d isoforms imported\n", nrow(aSwitchList$isoformFeatures)))

# ── 4. PreFilter ──────────────────────────────────────────────────
cat("\n[4/5] preFilter...\n")
aSwitchList <- preFilter(
  aSwitchList,
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.05,            # Less strict — catches more isoforms
  min.IF.prop = 0.25,         # Lower threshold to keep more isoforms
  removeSingleIsoformGenes = TRUE,
  alpha = 0.05,
  dIFcutoff = 0.1,
  quiet = FALSE
)

cat(sprintf("  → %d isoforms in %d genes after filter\n",
            nrow(aSwitchList$isoformFeatures),
            length(unique(aSwitchList$isoformFeatures$gene_id))))

# ── 5. DEXSeq test with parallel processing ──────────────────
cat("\n[5/5] isoformSwitchTestDEXSeq (all 5 comparisons, 60 cores)...\n")

# Override DEXSeq default BPPARAM to use parallel
# This modifies the function in the DEXSeq namespace
assignInNamespace("testForDEU", function(object, fullModel = design(object), 
    reducedModel = ~sample + exon, 
    BPPARAM = MulticoreParam(workers = 60), 
    fitType = c("DESeq2", "glmGamPoi")) {
    .Call(C_testForDEU, object, fullModel, reducedModel, BPPARAM, match.arg(fitType))
}, ns = "DEXSeq")

aSwitchList <- isoformSwitchTestDEXSeq(
  aSwitchList,
  alpha = 0.05,
  dIFcutoff = 0.1,
  reduceToSwitchingGenes = TRUE,
  showProgress = TRUE,
  quiet = FALSE
)

# ── Extract results per comparison ─────────────────────────────────
cat("\n════════════════════════════════════════════════════════════\n")
cat("  RESULTS\n")
cat("════════════════════════════════════════════════════════════\n\n")

all_features <- aSwitchList$isoformFeatures
results_list <- list()

for (i in 1:length(comp_names)) {
  cn <- comp_names[i]
  c1 <- comp_df$condition_1[i]
  c2 <- comp_df$condition_2[i]
  
  # Filter to this comparison
  comp_features <- all_features %>%
    filter(condition_1 == c1, condition_2 == c2)
  
  switches <- comp_features %>%
    filter(isoform_switch_q_value < 0.05, abs(dIF) >= 0.1)
  
  n_sig <- nrow(switches)
  n_genes <- length(unique(switches$gene_id))
  
  cat(sprintf("  %-25s: %4d switches (%3d genes)\n", cn, n_sig, n_genes))
  
  # Save per-comparison
  comp_out <- file.path(OUTDIR, cn)
  dir.create(comp_out, showWarnings = FALSE)
  write.csv(comp_features, file.path(comp_out, "all_isoform_features.csv"), row.names = FALSE)
  write.csv(switches, file.path(comp_out, "significant_switches.csv"), row.names = FALSE)
  
  if (n_sig > 0) {
    top <- switches %>% arrange(isoform_switch_q_value) %>% head(20)
    cat("\n    Top switches:\n")
    for (j in 1:nrow(top)) {
      cat(sprintf("    %2d. %-25s dIF=%+.4f q=%.2e iso=%s\n",
                  j, top$gene_id[j], top$dIF[j], top$isoform_switch_q_value[j],
                  top$isoform_id[j]))
    }
  }
  
  results_list[[cn]] <- list(
    comparison = cn,
    n_isoforms = nrow(comp_features),
    n_switches = n_sig,
    n_switch_genes = n_genes
  )
}

# Save full switch object
saveRDS(aSwitchList, file.path(OUTDIR, "switchAnalyzeRlist.rds"))

# ── Summary ────────────────────────────────────────────────────────
cat("\n  ────────\n")
cat("  SUMMARY\n")
cat("  ────────\n\n")
cat(sprintf("  %-30s %8s %8s\n", "Comparison", "Switches", "Genes"))
cat("  "), cat(rep("─", 48), sep=""); cat("\n")

overall <- data.frame()
for (cn in comp_names) {
  r <- results_list[[cn]]
  cat(sprintf("  %-30s %8d %8d\n", r$comparison, r$n_switches, r$n_switch_genes))
  overall <- rbind(overall, data.frame(
    Comparison = r$comparison,
    Switches = r$n_switches,
    Genes = r$n_switch_genes,
    stringsAsFactors = FALSE
  ))
}
write.csv(overall, file.path(OUTDIR, "summary.csv"), row.names = FALSE)

# Cross-comparison analysis
cat("\n\n  Cross-comparison: genes switching in multiple contrasts\n")
all_switch_genes <- list()
for (cn in comp_names) {
  sig_file <- file.path(OUTDIR, cn, "significant_switches.csv")
  if (file.exists(sig_file)) {
    df <- read.csv(sig_file)
    if (nrow(df) > 0) {
      all_switch_genes[[cn]] <- unique(df$gene_id)
    }
  }
}
all_genes <- unique(unlist(all_switch_genes))
gene_comp_count <- data.frame(
  gene = all_genes,
  n_comparisons = sapply(all_genes, function(g) sum(sapply(all_switch_genes, function(v) g %in% v))),
  comparisons = sapply(all_genes, function(g) 
    paste(names(all_switch_genes)[sapply(all_switch_genes, function(v) g %in% v)], collapse="; ")),
  stringsAsFactors = FALSE
)
gene_comp_count <- gene_comp_count[order(-gene_comp_count$n_comparisons), ]
cat(sprintf("  %d genes switch in ≥1 comparison\n", nrow(gene_comp_count)))
cat(sprintf("  %d genes switch in ≥2 comparisons\n", sum(gene_comp_count$n_comparisons >= 2)))
top_multi <- head(gene_comp_count[gene_comp_count$n_comparisons >= 2, ], 30)
if (nrow(top_multi) > 0) {
  write.csv(top_multi, file.path(OUTDIR, "multiswitch_genes.csv"), row.names = FALSE)
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("  ✅ ALL DONE\n")
cat(sprintf("  Results: %s/\n", OUTDIR))
cat("════════════════════════════════════════════════════════════\n")
