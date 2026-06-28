#!/usr/bin/env python3
"""
CTVT DESeq2 — Proper count matrix from prepDE.py
Gene-level: gene_count_matrix_prepDE.csv
Design: ~ dog_id + timepoint (9 dogs × 4 timepoints)
"""

import pandas as pd
import numpy as np
import os, warnings
warnings.filterwarnings('ignore')

# ======== Load prepDE count matrix ========
counts = pd.read_csv("gene_count_matrix_prepDE.csv", index_col=0)
counts = counts.round().astype(int)
print(f"Count matrix: {counts.shape[0]} genes × {counts.shape[1]} samples")

# ======== Filter ========
detected = (counts > 0).sum(axis=1)
keep = detected >= 9  # expressed in at least 9 samples (one full group)
count_filt = counts.loc[keep]
print(f"After filter (≥9 samples): {count_filt.shape[0]} genes")

# ======== Metadata ========
pheno = pd.read_csv("phenotype.csv", index_col=0)
pheno['dog_id'] = pheno['dog_id'].astype(str)
pheno['timepoint'] = pd.Categorical(pheno['timepoint'],
                                     categories=['Day0','Day2','Day6','Recovered'], ordered=True)
pheno = pheno.loc[count_filt.columns]
print(f"\nDesign: ~ dog_id + timepoint | {len(pheno)} samples")
print(pheno.groupby('timepoint').size().to_string())

# ======== DESeq2 ========
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=count_filt.T, metadata=pheno, design="~ dog_id + timepoint")
dds.deseq2()
print("DESeq2 fitted!")

# ======== 5 Contrasts ========
contrasts = [
    ("Day0_vs_Day2", "Day0", "Day2"),
    ("Day0_vs_Day6", "Day0", "Day6"),
    ("Day2_vs_Day6", "Day2", "Day6"),
    ("Day0_vs_Recovered", "Day0", "Recovered"),
    ("Day6_vs_Recovered", "Day6", "Recovered"),
]

outdir = "DE_results_final"
os.makedirs(outdir, exist_ok=True)

summary = []
for name, ref, alt in contrasts:
    stat = DeseqStats(dds, contrast=["timepoint", alt, ref], alpha=0.05)
    stat.summary()
    res = stat.results_df.dropna(subset=['padj']).sort_values('padj')
    
    sig = res[(res['padj'] < 0.05) & (abs(res['log2FoldChange']) >= 1)]
    up = (sig['log2FoldChange'] > 1).sum()
    dn = (sig['log2FoldChange'] < -1).sum()
    
    print(f"  {name:<25}: {len(sig):>6} DEGs ({up}↑ {dn}↓)")
    summary.append([name, len(sig), up, dn])

    res.to_csv(f"{outdir}/{name}_all.csv")
    if len(sig) > 0:
        sig.to_csv(f"{outdir}/{name}_significant.csv")
        top5 = sig.head(5)
        for i, (idx, row) in enumerate(top5.iterrows()):
            d = "↑" if row['log2FoldChange'] > 0 else "↓"
            print(f"    {i+1}. {idx}: {d} {2**row['log2FoldChange']:.1f}x (padj={row['padj']:.2e})")

# All sig genes
all_sig = set()
for name, _, _ in contrasts:
    try:
        df = pd.read_csv(f"{outdir}/{name}_significant.csv", index_col=0)
        all_sig.update(df.index)
    except: pass
print(f"\nUnique significant genes: {len(all_sig)}")

# Save final summary
summary_df = pd.DataFrame(summary, columns=["Contrast","DEGs","Up","Down"])
summary_df.to_csv(f"{outdir}/summary.csv", index=False)
count_filt.to_csv(f"{outdir}/filtered_counts.csv")
print(f"\nResults: {outdir}/")
print(summary_df.to_string(index=False))
