#!/usr/bin/env python3
"""
CTVT 36-sample: Coding potential analysis — FASTA-only version (no GTF parsing)
"""
import os, csv
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

BASE = "/home/vet/CTVT_raw_fastq_36samples"
OUT_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "coding_potential")
os.makedirs(OUT_DIR, exist_ok=True)

FA_NOVEL = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "novel_de_transcripts_only.fasta")

print("=" * 60)
print("  CTVT 36-sample — Coding Potential Analysis (FASTA-only)")
print("=" * 60)

# Load novel DE gene IDs
novel_de_genes = set()
with open(os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "novel_de_genes.txt")) as f:
    for line in f:
        novel_de_genes.add(line.strip())
print(f"  Novel DE genes: {len(novel_de_genes)}")

# ORF analysis functions
def find_longest_orf(seq, min_aa=30):
    seq = str(seq).upper().replace('N', '')
    if len(seq) < 3:
        return 0, 0.0, 0
    
    best_aa_len = 0
    best_nt_len = 0
    
    for strand, sequence in [(1, seq), (-1, str(Seq(seq).reverse_complement()))]:
        for frame in range(3):
            i = frame
            while i < len(sequence) - 2:
                codon = sequence[i:i+3]
                if codon == 'ATG':
                    aa_len = 0
                    j = i
                    while j < len(sequence) - 2:
                        c = sequence[j:j+3]
                        try:
                            aa = str(Seq(c).translate())
                            if aa == '*':
                                break
                        except:
                            break
                        aa_len += 1
                        j += 3
                    if aa_len > best_aa_len:
                        best_aa_len = aa_len
                        best_nt_len = j - i + 3
                i += 3
    
    coverage = best_nt_len / len(seq) if len(seq) > 0 else 0
    return best_aa_len, coverage, best_nt_len

def fickett_score(seq):
    seq = str(seq).upper()
    if len(seq) < 6:
        return 0.0
    n = len(seq)
    freqs = {'A': seq.count('A')/n, 'T': seq.count('T')/n, 
             'G': seq.count('G')/n, 'C': seq.count('C')/n}
    thirds = n // 3
    if thirds == 0:
        return 0.0
    pos = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for i, nt in enumerate(seq[:thirds*3]):
        if i % 3 == 0 and nt in pos:
            pos[nt] += 1
    for k in pos:
        pos[k] /= thirds
    score = sum(freqs[k] * pos[k] for k in 'ATGC') / 4
    return min(score, 1.0)

# Parse FASTA and analyze
print(f"\n  Analyzing coding potential of novel DE transcripts...")
results = []
seen = 0

for record in SeqIO.parse(FA_NOVEL, "fasta"):
    tx_id = record.id.split()[0]
    # Extract gene_id from transcript header (MSTRG.XXXX)
    # Transcripts are MSTRG.XXXXX.Y format, gene is MSTRG.XXXXX
    parts = tx_id.split('.')
    gene_id = '.'.join(parts[:-1]) if len(parts) > 1 else tx_id
    
    seq = str(record.seq)
    seq_len = len(seq)
    if seq_len < 30:
        continue
    
    orf_aa, orf_cov, orf_nt = find_longest_orf(seq)
    fickett = fickett_score(seq)
    
    # Classification
    if orf_aa == 0 or orf_aa < 30:
        cls = "high_confidence_lncRNA"
    elif orf_aa < 100 and orf_cov < 0.1 and fickett < 0.1:
        cls = "high_confidence_lncRNA"
    elif orf_aa < 100:
        cls = "likely_lncRNA"
    elif orf_aa < 300 and orf_cov < 0.3:
        cls = "likely_lncRNA"
    elif orf_aa < 300:
        cls = "ambiguous"
    else:
        cls = "potential_coding"
    
    results.append({
        'transcript_id': tx_id,
        'gene_id': gene_id,
        'seq_length': seq_len,
        'longest_ORF_aa': orf_aa,
        'ORF_coverage': round(orf_cov, 3),
        'ORF_nt_length': orf_nt,
        'fickett_score': round(fickett, 4),
        'classification': cls
    })
    
    seen += 1
    if seen % 2000 == 0:
        print(f"    Processed {seen} transcripts...")

print(f"\n  Total novel DE transcripts analyzed: {len(results)}")

df = pd.DataFrame(results)
df.to_csv(os.path.join(OUT_DIR, "coding_potential_results.tsv"), sep='\t', index=False)

for cls_name in ['high_confidence_lncRNA', 'likely_lncRNA', 'ambiguous', 'potential_coding']:
    sub = df[df['classification'] == cls_name]
    sub.to_csv(os.path.join(OUT_DIR, f"{cls_name}.tsv"), sep='\t', index=False)
    print(f"  {cls_name}: {len(sub)} transcripts")

lncrna_genes = set(df[df['classification'].isin(['high_confidence_lncRNA', 'likely_lncRNA'])]['gene_id'])
with open(os.path.join(OUT_DIR, "lncRNA_candidate_gene_ids.txt"), 'w') as f:
    for g in sorted(lncrna_genes):
        f.write(g + '\n')

print(f"\n  Combined lncRNA candidates (HC + Likely): {len(df[df['classification'].isin(['high_confidence_lncRNA', 'likely_lncRNA'])])} transcripts")
print(f"  Unique lncRNA gene IDs: {len(lncrna_genes)}")
print("\nStep 2 complete!")