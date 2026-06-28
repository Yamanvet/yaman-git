#!/usr/bin/env Rscript
# ====================================================================
# CTVT: Enrichment analysis — 36-sample protein-coding DEGs
# KEGG, GO (BP/CC/MF), Reactome, MSigDB Hallmark
# Using clusterProfiler + enrichR
# ====================================================================

Sys.setenv(R_MAX_VSIZE = "250e9")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(enrichR)
})

BASE <- "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/protein_coding"
OUTDIR <- "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/enrichment"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Load MSTRG → symbol mapping
mstrg_map <- read.csv("/home/vet/CTVT_raw_fastq_36samples/mstrg_to_symbol.csv", stringsAsFactors = FALSE)
symbol_map <- setNames(mstrg_map$symbol, mstrg_map$gene)

contrasts <- c("Day2_vs_Day0","Day6_vs_Day0","Day6_vs_Day2","Recovered_vs_Day0","Recovered_vs_Day6")

# EnrichR databases
enrichr_dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", 
                  "GO_Molecular_Function_2023", "KEGG_2021_Human", 
                  "Reactome_2022", "MSigDB_Hallmark_2020")

for (cname in contrasts) {
  cat(sprintf("\n%s %s %s\n", strrep("═", 20), cname, strrep("═", 20)))
  
  f <- file.path(BASE, paste0(cname, "_DEGs.csv"))
  if (!file.exists(f)) {
    cat("  ⚠️ File not found, skipping\n")
    next
  }
  
  degs <- read.csv(f, stringsAsFactors = FALSE)
  # Map MSTRG to symbols
  degs$symbol <- symbol_map[degs$gene]
  degs <- degs[!is.na(degs$symbol) & degs$symbol != "", ]
  
  # Split up/down
  degs_up <- degs[degs$log2FoldChange > 0, ]
  degs_down <- degs[degs$log2FoldChange < 0, ]
  degs_all <- degs
  
  cat(sprintf("  Total DEGs: %d (up=%d, down=%d)\n", nrow(degs_all), nrow(degs_up), nrow(degs_down)))
  
  # Convert symbols to Entrez IDs for clusterProfiler
  all_symbols <- degs_all$symbol
  up_symbols <- degs_up$symbol
  down_symbols <- degs_down$symbol
  
  all_entrez <- bitr(all_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
  up_entrez <- bitr(up_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  down_entrez <- bitr(down_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
  
  cat(sprintf("  Mapped: all=%d, up=%d, down=%d\n", 
    nrow(all_entrez), nrow(up_entrez), nrow(down_entrez)))
  
  cdir <- file.path(OUTDIR, cname)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  
  # ── clusterProfiler: KEGG ──────────────────────────────────────
  for (direction in c("up", "down", "all")) {
    genes <- switch(direction,
      up = up_entrez$ENTREZID,
      down = down_entrez$ENTREZID,
      all = all_entrez$ENTREZID
    )
    
    if (length(genes) < 5) {
      cat(sprintf("  KEGG %s: too few genes (%d), skipping\n", direction, length(genes)))
      next
    }
    
    kk <- tryCatch({
      enrichKEGG(gene = genes, organism = "hsa", pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, minGSSize = 5)
    }, error = function(e) NULL)
    
    if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
      write.csv(as.data.frame(kk), file.path(cdir, paste0("KEGG_", direction, ".csv")), row.names = FALSE)
      cat(sprintf("  KEGG %s: %d pathways\n", direction, nrow(as.data.frame(kk))))
    } else {
      cat(sprintf("  KEGG %s: 0 pathways\n", direction))
    }
  }
  
  # ── clusterProfiler: GO ────────────────────────────────────────
  for (ont in c("BP", "CC", "MF")) {
    for (direction in c("up", "down", "all")) {
      genes <- switch(direction,
        up = up_entrez$ENTREZID,
        down = down_entrez$ENTREZID,
        all = all_entrez$ENTREZID
      )
      
      if (length(genes) < 5) next
      
      ego <- tryCatch({
        enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = ont,
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 5,
                  readable = TRUE)
      }, error = function(e) NULL)
      
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), file.path(cdir, paste0("GO_", ont, "_", direction, ".csv")), row.names = FALSE)
        cat(sprintf("  GO %s %s: %d terms\n", ont, direction, nrow(as.data.frame(ego))))
      } else {
        cat(sprintf("  GO %s %s: 0 terms\n", ont, direction))
      }
    }
  }
  
  # ── EnrichR: Reactome, MSigDB Hallmark ─────────────────────────
  for (direction in c("up", "down", "all")) {
    genes <- switch(direction,
      up = up_symbols,
      down = down_symbols,
      all = all_symbols
    )
    
    if (length(genes) < 5) next
    
    enr <- tryCatch({
      enrichr(genes, databases = enrichr_dbs)
    }, error = function(e) NULL)
    
    if (!is.null(enr)) {
      for (db in names(enr)) {
        df <- enr[[db]]
        df <- df[df$Adjusted.P.value < 0.05, ]
        db_clean <- gsub("[^a-zA-Z0-9]", "_", db)
        if (nrow(df) > 0) {
          write.csv(df, file.path(cdir, paste0(db_clean, "_", direction, ".csv")), row.names = FALSE)
          cat(sprintf("  EnrichR %s %s: %d terms\n", db, direction, nrow(df)))
        }
      }
    }
  }
}

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("  36-sample enrichment complete!\n")
cat("═══════════════════════════════════════════════════════════════\n")
