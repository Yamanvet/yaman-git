#!/usr/bin/env python3
"""
CTVT: Complete 36-sample downstream analysis
Mirrors the 33-sample pipeline: gseapy enrichment (full reports) + volcano plots + heatmaps + GSEA preranked
"""

import os, csv, numpy as np, pandas as pd
import gseapy as gp
from gseapy.plot import barplot, dotplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

BASE = "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36"
OUTDIR = "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/figures"
os.makedirs(OUTDIR, exist_ok=True)

# Load gene name mapping
mstrg_map = {}
with open("/home/vet/CTVT_raw_fastq_36samples/mstrg_to_gene_name.csv") as f:
    for row in csv.DictReader(f):
        mstrg_map[row['MSTRG_ID']] = row['gene_name']

contrasts = ["Day2_vs_Day0", "Day6_vs_Day0", "Day6_vs_Day2", "Recovered_vs_Day0", "Recovered_vs_Day6"]

# ── 1. Load all DEG data ──────────────────────────
print("Loading DEG data...")
all_degs = {}
for cname in contrasts:
    f = os.path.join(BASE, "protein_coding", f"{cname}_DEGs.csv")
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, index_col=0)
    df['symbol'] = df['gene'].map(mstrg_map).fillna(df['gene'])
    df = df.drop_duplicates(subset='symbol')
    all_degs[cname] = df

# ── 2. Volcano plots ──────────────────────────
print("Generating volcano plots...")
fig, axes = plt.subplots(2, 3, figsize=(22, 14))
axes = axes.flatten()

for i, cname in enumerate(contrasts):
    ax = axes[i]
    df = all_degs[cname]
    
    sig_up = (df['padj'] < 0.05) & (df['log2FoldChange'] > 0)
    sig_down = (df['padj'] < 0.05) & (df['log2FoldChange'] < 0)
    nonsig = ~(sig_up | sig_down)
    
    ax.scatter(df.loc[nonsig, 'log2FoldChange'], -np.log10(df.loc[nonsig, 'padj']),
               c='gray', s=5, alpha=0.3, label=f'NS ({nonsig.sum()})')
    ax.scatter(df.loc[sig_up, 'log2FoldChange'], -np.log10(df.loc[sig_up, 'padj']),
               c='red', s=10, alpha=0.6, label=f'UP ({sig_up.sum()})')
    ax.scatter(df.loc[sig_down, 'log2FoldChange'], -np.log10(df.loc[sig_down, 'padj']),
               c='blue', s=10, alpha=0.6, label=f'DOWN ({sig_down.sum()})')
    
    # Label top genes
    top_up = df[sig_up].nsmallest(5, 'padj')
    top_down = df[sig_down].nsmallest(5, 'padj')
    for _, row in pd.concat([top_up, top_down]).iterrows():
        ax.annotate(row['symbol'], (row['log2FoldChange'], -np.log10(row['padj'])),
                   fontsize=6, ha='center', va='bottom')
    
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.5)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax.set_title(cname, fontsize=11, fontweight='bold')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(padj)')
    ax.legend(fontsize=7, loc='upper right')

# Hide empty subplot
axes[5].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'volcano_plots_36sample.png'), dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved volcano_plots_36sample.png")

# ── 3. Heatmap of top 50 DEGs (Recovered vs Day0) ──────────────────────────
print("Generating heatmap...")

# Load normalized counts
norm_file = os.path.join(BASE, "protein_coding", "normalized_counts.csv")
if os.path.exists(norm_file):
    norm_counts = pd.read_csv(norm_file, index_col=0)
    # Map gene IDs to symbols
    new_idx = [mstrg_map.get(g, g) for g in norm_counts.index]
    norm_counts.index = new_idx
    
    # Use the biggest contrast for the heatmap
    cname = "Recovered_vs_Day0"
    deg_df = all_degs[cname]
    top50 = deg_df.nsmallest(50, 'padj')
    top50_genes = top50['symbol'].tolist()
    
    # Filter to genes present in normalized counts
    top50_genes = [g for g in top50_genes if g in norm_counts.index]
    
    if len(top50_genes) >= 10:
        heat_data = norm_counts.loc[top50_genes]
        heat_data = heat_data.dropna()
        
        # Z-score normalize
        heat_z = heat_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        
        # Sort columns by timepoint
        sample_order = []
        for tp in ['Day0', 'Day2', 'Day6', 'Recovered']:
            samples = [s for s in heat_z.columns if tp in s or (tp == 'Day0' and 'D0' in s)]
            sample_order.extend(sorted(samples))
        # Fallback: use all columns
        if not sample_order:
            sample_order = list(heat_z.columns)
        
        sample_order = [s for s in sample_order if s in heat_z.columns]
        if len(sample_order) < 2:
            sample_order = list(heat_z.columns)
        
        heat_z = heat_z[sample_order]
        
        fig, ax = plt.subplots(figsize=(16, 20))
        sns.heatmap(heat_z, cmap='RdBu_r', center=0, ax=ax,
                   xticklabels=True, yticklabels=True,
                   linewidths=0.5, linecolor='white')
        ax.set_title(f'Top {len(top50_genes)} DEGs: {cname} (36 samples)', fontsize=14)
        ax.set_xlabel('Sample')
        ax.set_ylabel('Gene')
        plt.tight_layout()
        plt.savefig(os.path.join(OUTDIR, 'heatmap_top50_36sample.png'), dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  Saved heatmap_top50_36sample.png")
    else:
        print("  Not enough genes for heatmap")
else:
    print(f"  Normalized counts file not found: {norm_file}")

# ── 4. gseapy enrichment (full reports, all terms) ──────────────────────────
print("\nRunning gseapy enrichment (full reports)...")

enrichr_dbs = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    "BioCarta_2016",
    "Elsevier_Pathway_Collection"
]

for cname in contrasts:
    if cname not in all_degs:
        continue
    df = all_degs[cname]
    
    for direction in ["up", "down", "all"]:
        if direction == "up":
            genes = df[df['log2FoldChange'] > 0]['symbol'].tolist()
        elif direction == "down":
            genes = df[df['log2FoldChange'] < 0]['symbol'].tolist()
        else:
            genes = df['symbol'].tolist()
        
        if len(genes) < 5:
            print(f"  {cname} {direction}: too few genes ({len(genes)}), skip")
            continue
        
        outdir = os.path.join(BASE, "enrichment_gseapy", f"{cname}_{direction}")
        os.makedirs(outdir, exist_ok=True)
        
        print(f"  {cname} {direction} ({len(genes)} genes)...")
        try:
            result = gp.enrichr(gene_list=genes, gene_sets=enrichr_dbs,
                               organism='human', outdir=outdir, no_plot=True,
                               cutoff=1.0)  # Get ALL terms, not just sig
            print(f"    Done")
        except Exception as e:
            print(f"    Error: {e}")

# ── 5. GSEA preranked (shrunken LFC) ──────────────────────────
print("\nRunning GSEA preranked...")

for cname in contrasts:
    if cname not in all_degs:
        continue
    df = all_degs[cname]
    
    # Create ranked gene list (shrunken LFC)
    rnk = df[['symbol', 'shrunkenLog2FC']].copy()
    rnk = rnk.dropna()
    rnk = rnk.drop_duplicates(subset='symbol')
    rnk = rnk.sort_values('shrunkenLog2FC', ascending=False)
    rnk.columns = ['gene', 'score']
    
    outdir = os.path.join(BASE, "gsea_preranked", cname)
    os.makedirs(outdir, exist_ok=True)
    
    # Save rnk file
    rnk.to_csv(os.path.join(outdir, "ranked_genes.rnk"), sep='\t', index=False, header=False)
    
    print(f"  GSEA {cname} ({len(rnk)} genes)...")
    try:
        gs = gp.prerank(rnk=rnk, gene_sets='KEGG_2021_Human',
                       organism='human', outdir=outdir,
                       permutation_num=1000, seed=42, no_plot=True,
                       min_size=5, max_size=500)
        print(f"    KEGG done")
    except Exception as e:
        print(f"    KEGG error: {e}")
    
    # Also run Hallmark
    try:
        gs2 = gp.prerank(rnk=rnk, gene_sets='MSigDB_Hallmark_2020',
                        organism='human', outdir=os.path.join(outdir, "hallmark"),
                        permutation_num=1000, seed=42, no_plot=True,
                        min_size=5, max_size=500)
        print(f"    Hallmark done")
    except Exception as e:
        print(f"    Hallmark error: {e}")

print("\nAll analyses complete!")