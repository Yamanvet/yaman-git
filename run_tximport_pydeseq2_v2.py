#!/usr/bin/env python3
"""
CTVT: tximport + PyDESeq2 — Gene-level DE analysis v2
Drops 4 outlier samples: T1A_S1, T1B_S2, T8C_S10, T9N_S26
Runs PyDESeq2 with ~ dog_id + timepoint design (n=8 dogs, 32 samples)
"""

import os, warnings, subprocess, tempfile
import pandas as pd
import numpy as np

warnings.filterwarnings('ignore')

BASE = "/home/vet/CTVT_raw_fastq_36samples"
QUANT_DIR = "/mnt/data8tb/salmon_quant_final"
OUTDIR = os.path.join(BASE, "deseq2_v4_tximport_outlier_removed_v3")
GTF_PATH = os.path.join(BASE, "stringtie_v4", "merged_stringtie_v4.gtf")
os.makedirs(OUTDIR, exist_ok=True)

# Samples to drop (keeping T1A_S1 this run)
DROP_SAMPLES = {"T1B_S2", "T8C_S10", "T9N_S26"}

print("="*60)
print("  CTVT: tximport + PyDESeq2 — Outlier-Removed v2")
print("  Dropped: T1B_S2, T8C_S10, T9N_S26 (kept T1A_S1)")
print("="*60)

# ── 1. Phenotype ─────────────────────────────────────
samples = sorted([s for s in os.listdir(QUANT_DIR) if s not in DROP_SAMPLES])
pheno = pd.DataFrame({"sample_id": samples})
pheno["dog_id"] = pheno["sample_id"].str.extract(r"(T\d+)")
timepoint_map = {"A":"Day0","B":"Day2","C":"Day6","N":"Recovered"}
pheno["timepoint"] = pheno["sample_id"].str[2].map(timepoint_map)
pheno.to_csv(os.path.join(OUTDIR, "phenotype.csv"), index=False)
print(f"\nDesign: {len(pheno)} samples (after dropping 4)")
print(pheno.groupby("timepoint").size().to_string())

# ── 2. tximport via R with rtracklayer ───────────────
print("\n[1/4] Running tximport (Salmon → gene counts)...")

r_script = f'''
suppressPackageStartupMessages({{
  library(tximport, lib.loc="~/R/library")
  library(rtracklayer, lib.loc="~/R/library")
}})

quant_dir <- "{QUANT_DIR}"
samples <- list.files(quant_dir)
files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples

gtf <- import("{GTF_PATH}")
tx_data <- as.data.frame(gtf[gtf$type == "transcript"])
tx2gene <- unique(tx_data[, c("transcript_id", "gene_id")])
cat(nrow(tx2gene), "transcript-gene pairs loaded\\n")

txi <- tximport(files, type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM", dropInfReps=TRUE)
write.csv(txi$counts, "{OUTDIR}/gene_count_matrix_full.csv")
write.csv(txi$abundance, "{OUTDIR}/gene_tpm_matrix_full.csv")
cat("Done:", nrow(txi$counts), "genes x", ncol(txi$counts), "samples\\n")
'''

with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
    f.write(r_script)
    rfile = f.name

subprocess.run(["Rscript", rfile], check=True)
os.unlink(rfile)

# ── 3. Load, drop outliers, filter ───────────────────
print("\n[2/4] Loading, dropping outliers, and filtering...")
counts_full = pd.read_csv(f"{OUTDIR}/gene_count_matrix_full.csv", index_col=0)
# Drop outlier samples
keep_cols = [c for c in counts_full.columns if c not in DROP_SAMPLES]
counts = counts_full[keep_cols].round().astype(int)
print(f"  After dropping outliers: {counts.shape[1]} samples")

detected = (counts > 0).sum(axis=1)
# Filter: detected in >= min_group_size (now 8 per group, minus T1/T8/T9 losses)
# Day0: 8, Day2: 8, Day6: 8, Recovered: 8
counts_filt = counts.loc[detected >= 8]
print(f"  Raw: {counts.shape[0]} genes × {counts.shape[1]} samples")
print(f"  After filter (≥8 samples): {counts_filt.shape[0]} genes")

# ── 4. PyDESeq2 ──────────────────────────────────────
print("\n[3/4] Running PyDESeq2...")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

pheno_filt = pheno.set_index("sample_id").loc[counts_filt.columns]
pheno_filt["dog_id"] = pheno_filt["dog_id"].astype(str)
pheno_filt["timepoint"] = pd.Categorical(pheno_filt["timepoint"],
    categories=["Day0","Day2","Day6","Recovered"], ordered=True)

dds = DeseqDataSet(counts=counts_filt.T, metadata=pheno_filt,
                   design="~ dog_id + timepoint", n_cpus=16)
dds.deseq2()
print("  ✅ DESeq2 fitted!")

# ── 5. All 5 comparisons ─────────────────────────────
print("\n[4/4] Running 5 pairwise comparisons...")
comparisons = [
    ("Day2_vs_Day0", "Day2", "Day0"),
    ("Day6_vs_Day0", "Day6", "Day0"),
    ("Recovered_vs_Day0", "Recovered", "Day0"),
    ("Day6_vs_Day2", "Day6", "Day2"),
    ("Recovered_vs_Day6", "Recovered", "Day6"),
]

results = []
for name, c1, c2 in comparisons:
    print(f"\n  {name} ({c1} vs {c2})...")
    stat = DeseqStats(dds, contrast=["timepoint", c1, c2], n_cpus=16)
    stat.summary()
    
    all_res = stat.results_df.copy()
    sig = all_res[all_res.padj < 0.05].copy()
    
    all_res.to_csv(f"{OUTDIR}/{name}_all.csv")
    sig.to_csv(f"{OUTDIR}/{name}_significant.csv")
    
    n_sig = len(sig)
    n_up = sum(sig.log2FoldChange > 0) if n_sig > 0 else 0
    n_down = sum(sig.log2FoldChange < 0) if n_sig > 0 else 0
    results.append({"Comparison":name, "DE_genes":n_sig, "Up":n_up, "Down":n_down})
    print(f"    {n_sig} DE genes ({n_up} up, {n_down} down)")

# ── Summary ──────────────────────────────────────────
print("\n" + "="*60)
print("  RESULTS SUMMARY (Outlier-Removed)")
print("="*60)
summary = pd.DataFrame(results)
print(summary.to_string(index=False))
summary.to_csv(f"{OUTDIR}/DE_summary.csv", index=False)

# Compare with original
orig_summary = pd.read_csv(os.path.join(BASE, "deseq2_v4_tximport", "DE_summary.csv"))
print("\n  Comparison with original (all 36 samples):")
for _, row in orig_summary.iterrows():
    new_row = summary[summary.Comparison == row.Comparison]
    if len(new_row) > 0:
        diff = int(new_row.DE_genes.values[0]) - int(row.DE_genes)
        print(f"    {row.Comparison}: {row.DE_genes} → {new_row.DE_genes.values[0]} ({'+' if diff >= 0 else ''}{diff})")

print(f"\nAll results: {OUTDIR}/")