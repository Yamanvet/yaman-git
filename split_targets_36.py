#!/usr/bin/env python3
"""
CTVT 36-sample: Split DEGs into Known (with gene names) and Novel (lncRNA candidates)
Also generates ranked gene lists for GSEA and pathway summaries
"""

import os, csv, re
import pandas as pd
import numpy as np

BASE = "/home/vet/CTVT_raw_fastq_36samples"
DEG_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "protein_coding")
OUT_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "known_genes_pathways")
NOVEL_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "novel_lncRNA_candidates")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(NOVEL_DIR, exist_ok=True)

# Load MSTRG → gene_name mapping from merged GTF
print("Building gene annotation map from merged GTF...")
mstrg_to_name = {}
mstrg_to_ref = {}
gene_type = {}

GTF_PATH = os.path.join(BASE, "stringtie_v4", "merged_stringtie_v4.gtf")
with open(GTF_PATH) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        attr = fields[8]
        
        gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
        if not gene_id_match:
            continue
        gene_id = gene_id_match.group(1)
        
        gene_name_match = re.search(r'gene_name "([^"]+)"', attr)
        ref_gene_match = re.search(r'ref_gene_id "([^"]+)"', attr)
        
        gene_name = gene_name_match.group(1) if gene_name_match else None
        ref_gene = ref_gene_match.group(1) if ref_gene_match else None
        
        if ref_gene or gene_name:
            mstrg_to_name[gene_id] = gene_name or ref_gene
            mstrg_to_ref[gene_id] = ref_gene
            gene_type[gene_id] = "known"
        else:
            gene_type[gene_id] = "novel"

n_known = sum(1 for v in gene_type.values() if v == "known")
n_novel = sum(1 for v in gene_type.values() if v == "novel")
print(f"  {n_known} MSTRG genes with known ref names")
print(f"  {n_novel} MSTRG genes that are novel (no ref_gene_id)")

contrasts = ["Day2_vs_Day0", "Day6_vs_Day0", "Day6_vs_Day2", "Recovered_vs_Day0", "Recovered_vs_Day6"]

for cname in contrasts:
    deg_file = os.path.join(DEG_DIR, f"{cname}_DEGs.csv")
    all_file = os.path.join(DEG_DIR, f"{cname}_all_results.csv")
    
    if not os.path.exists(deg_file):
        print(f"  Skipping {cname} — DEG file not found")
        continue
    
    # Load DEGs
    deg_df = pd.read_csv(deg_file, index_col=0)
    
    # Add gene_name, ref_gene_id, gene_type
    deg_df['gene_name'] = deg_df['gene'].map(mstrg_to_name).fillna(deg_df['gene'])
    deg_df['ref_gene_id'] = deg_df['gene'].map(mstrg_to_ref)
    deg_df['gene_type'] = deg_df['gene'].map(gene_type).fillna('unknown')
    
    # Split known vs novel
    known_df = deg_df[deg_df['gene_type'] == 'known'].copy()
    novel_df = deg_df[deg_df['gene_type'] == 'novel'].copy()
    
    # Save
    known_df.to_csv(os.path.join(OUT_DIR, f"{cname}_known_genes.csv"))
    novel_df.to_csv(os.path.join(NOVEL_DIR, f"{cname}_novel_MSTRG.csv"))
    
    # Also save ALL results split
    all_df = pd.read_csv(all_file, index_col=0)
    all_df['gene_name'] = all_df['gene'].map(mstrg_to_name).fillna(all_df['gene'])
    all_df['ref_gene_id'] = all_df['gene'].map(mstrg_to_ref)
    all_df['gene_type'] = all_df['gene'].map(gene_type).fillna('unknown')
    
    all_known = all_df[all_df['gene_type'] == 'known'].copy()
    all_novel = all_df[all_df['gene_type'] == 'novel'].copy()
    
    all_known.to_csv(os.path.join(OUT_DIR, f"{cname}_all_known.csv"))
    all_novel.to_csv(os.path.join(NOVEL_DIR, f"{cname}_all_novel.csv"))
    
    # Create ranked gene list (for GSEA) — known genes only
    known_all = all_known[['gene_name', 'shrunkenLog2FC']].dropna(subset=['gene_name'])
    known_all = known_all.drop_duplicates(subset='gene_name', keep='first')
    known_all = known_all.sort_values('shrunkenLog2FC', ascending=False)
    known_all.columns = ['gene', 'score']
    known_all.to_csv(os.path.join(OUT_DIR, f"{cname}_ranked_known.rnk"), sep='\t', index=False, header=False)
    
    print(f"  {cname}: {len(known_df)} known DEGs, {len(novel_df)} novel DEGs, {len(all_known)} known total, {len(all_novel)} novel total")

print("\nDone! Known genes and novel lncRNA candidates saved.")