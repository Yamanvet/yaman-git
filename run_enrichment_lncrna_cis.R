#!/usr/bin/env Rscript
# CTVT: Enrichment on lncRNA cis-target genes (36 samples)

Sys.setenv(R_MAX_VSIZE = "250e9")
suppressPackageStartupMessages({
  library(enrichR)
})

setEnrichrSite("Enrichr")

BASE <- "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/lncrna_known"
OUTDIR <- "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/enrichment_lncrna"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

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
  cat(sprintf("\n%s %s %s\n", strrep("=", 25), cname, strrep("=", 25)))
  
  cis_file <- file.path(BASE, paste0(cname, "_lncrna_cis_targets.csv"))
  if (!file.exists(cis_file)) {
    cat("  No cis-target file, skipping\n")
    next
  }
  
  cis <- read.csv(cis_file, stringsAsFactors = FALSE)
  genes <- unique(cis$cis_gene_name)
  genes <- genes[!is.na(genes) & genes != ""]
  
  cat(sprintf("  Cis-target genes: %d\n", length(genes)))
  
  if (length(genes) < 5) {
    cat("  Too few genes, skipping\n")
    next
  }
  
  cdir <- file.path(OUTDIR, paste0(cname, "_cis_targets"))
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("  Running EnrichR on %d cis-target genes...\n", length(genes)))
  
  enr <- tryCatch({
    enrichr(genes, databases = enrichr_dbs)
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
      write.csv(df, file.path(cdir, paste0(db_clean, "_cis_targets.csv")), row.names = FALSE)
      cat(sprintf("    %s: %d terms\n", db, nrow(df)))
    } else {
      cat(sprintf("    %s: 0 terms\n", db))
    }
  }
}

cat("\n  lncRNA cis-target enrichment complete!\n")