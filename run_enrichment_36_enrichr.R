#!/usr/bin/env Rscript
# ====================================================================
# CTVT: Enrichment analysis — 36-sample protein-coding DEGs
# EnrichR only (no clusterProfiler needed)
# KEGG, GO (BP/CC/MF), Reactome, MSigDB Hallmark, BioCarta, Elsevier
# ====================================================================

Sys.setenv(R_MAX_VSIZE = "250e9")
suppressPackageStartupMessages({
  library(enrichR)
  library(readr)
})

# Force live access
enrichr:::.onAttach()

BASE <- "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/protein_coding"
OUTDIR <- "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/enrichment"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Load MSTRG → symbol mapping
mstrg_map <- read.csv("/home/vet/CTVT_raw_fastq_36samples/mstrg_to_symbol.csv", stringsAsFactors = FALSE)
symbol_map <- setNames(mstrg_map$symbol, mstrg_map$gene)

contrasts <- c("Day2_vs_Day0","Day6_vs_Day0","Day6_vs_Day2","Recovered_vs_Day0","Recovered_vs_Day6")

enrichr_dbs <- c(
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023", 
  "GO_Molecular_Function_2023",
  "KEGG_2021_Human",
  "Reactome_2022",
  "MSigDB_Hallmark_2020",
  "BioCarta_2019",
  "Elsevier_Pathway_Collection"
)

for (cname in contrasts) {
  cat(sprintf("\n%s %s %s\n", strrep("═", 20), cname, strrep("═", 20)))
  
  f <- file.path(BASE, paste0(cname, "_DEGs.csv"))
  if (!file.exists(f)) {
    cat("  ⚠️ File not found, skipping\n")
    next
  }
  
  degs <- read.csv(f, stringsAsFactors = FALSE)
  degs$symbol <- symbol_map[degs$gene]
  degs <- degs[!is.na(degs$symbol) & degs$symbol != "", ]
  
  # Remove duplicates (keep lowest padj)
  degs <- degs[order(degs$padj), ]
  degs <- degs[!duplicated(degs$symbol), ]
  
  degs_up <- degs[degs$log2FoldChange > 0, ]
  degs_down <- degs[degs$log2FoldChange < 0, ]
  
  cat(sprintf("  Total DEGs: %d (up=%d, down=%d)\n", nrow(degs), nrow(degs_up), nrow(degs_down)))
  
  cdir <- file.path(OUTDIR, cname)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  
  for (direction in c("up", "down", "all")) {
    genes <- switch(direction,
      up = degs_up$symbol,
      down = degs_down$symbol,
      all = degs$symbol
    )
    
    if (length(genes) < 5) {
      cat(sprintf("  %s %s: too few genes (%d), skipping\n", direction, "(?)", length(genes)))
      next
    }
    
    cat(sprintf("  Running EnrichR for %s (%d genes)...\n", direction, length(genes)))
    
    enr <- tryCatch({
      enrichr(genes, databases = enrichr_dbs)
    }, error = function(e) {
      cat(sprintf("  ⚠️ EnrichR error: %s\n", e$message))
      NULL
    })
    
    if (is.null(enr)) next
    
    for (db in names(enr)) {
      df <- enr[[db]]
      df <- df[df$Adjusted.P.value < 0.05, ]
      db_clean <- gsub("[^a-zA-Z0-9]", "_", db)
      
      if (nrow(df) > 0) {
        write.csv(df, file.path(cdir, paste0(db_clean, "_", direction, ".csv")), row.names = FALSE)
        cat(sprintf("    %s %s: %d terms\n", db, direction, nrow(df)))
      } else {
        cat(sprintf("    %s %s: 0 terms\n", db, direction))
      }
    }
  }
}

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  36-sample enrichment complete!\n")
cat("═══════════════════════════════════════════════════════════════\n")