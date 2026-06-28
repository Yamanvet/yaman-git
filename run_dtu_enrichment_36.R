#!/usr/bin/env Rscript
# CTVT: DTU gene-level enrichment (36 samples)

Sys.setenv(R_MAX_VSIZE = "250e9")
suppressPackageStartupMessages({
  library(enrichR)
})

setEnrichrSite("Enrichr")

BASE <- "/home/vet/CTVT_raw_fastq_36samples/diffsplice_v7_all36"
OUTDIR <- "/home/vet/CTVT_raw_fastq_36samples/diffsplice_v7_all36/enrichment"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# MSTRG → gene_name mapping
mstrg_map <- read.csv("/home/vet/CTVT_raw_fastq_36samples/mstrg_to_gene_name.csv", stringsAsFactors = FALSE)
symbol_map <- setNames(mstrg_map$gene_name, mstrg_map$MSTRG_ID)

contrasts <- c("Day0_vs_Day2","Day0_vs_Day6","Day0_vs_Recovered","Day2_vs_Day6","Day6_vs_Recovered")

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
  cat(sprintf("\n%s %s %s\n", strrep("=", 25), cname, strrep("=", 25)))
  
  gene_file <- file.path(BASE, cname, "significant_genes_q005.csv")
  if (!file.exists(gene_file)) { cat("  No gene file\n"); next }
  
  genes_df <- read.csv(gene_file, stringsAsFactors = FALSE, row.names = 1)
  genes_df$symbol <- symbol_map[genes_df$GeneID]
  genes_df$symbol[is.na(genes_df$symbol)] <- genes_df$GeneID[is.na(genes_df$symbol)]
  genes <- unique(genes_df$symbol[genes_df$symbol != ""])
  
  # Also filter to protein-coding only
  biotype_map <- read.csv("/home/vet/CTVT_raw_fastq_36samples/tx2gene_protein_coding.csv", stringsAsFactors = FALSE)
  pc_genes <- unique(biotype_map$GENEID)
  pc_genes_names <- unique(symbol_map[pc_genes[pc_genes %in% names(symbol_map)]])
  
  genes_pc <- genes[genes %in% pc_genes_names]
  
  cat(sprintf("  Total genes: %d, protein-coding: %d\n", length(genes), length(genes_pc)))
  
  if (length(genes_pc) < 5) {
    cat("  Too few PC genes for enrichment\n")
    next
  }
  
  cdir <- file.path(OUTDIR, cname)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("  EnrichR on %d protein-coding DTU genes...\n", length(genes_pc)))
  
  enr <- tryCatch({
    enrichr(genes_pc, databases = enrichr_dbs)
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    NULL
  })
  
  if (is.null(enr)) next
  
  for (db in names(enr)) {
    df <- enr[[db]]
    df <- df[df$Adjusted.P.value < 0.05, ]
    db_clean <- gsub("[^a-zA-Z0-9]", "_", db)
    
    if (nrow(df) > 0) {
      write.csv(df, file.path(cdir, paste0(db_clean, "_dtu_genes.csv")), row.names = FALSE)
      cat(sprintf("    %s: %d terms\n", db, nrow(df)))
    } else {
      cat(sprintf("    %s: 0 terms\n", db))
    }
  }
  
  rm(enr); gc()
}

cat("\n  DTU enrichment complete!\n")