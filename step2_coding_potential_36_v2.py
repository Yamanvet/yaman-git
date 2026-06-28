#!/usr/bin/env python3
"""
CTVT 36-sample: Step 2 — Coding potential for novel DE lncRNAs only
Optimized: filter FASTA to novel DE transcripts first, then analyze
"""
import os, re, csv
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd

BASE = "/home/vet/CTVT_raw_fastq_36samples"
OUT_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "coding_potential")
os.makedirs(OUT_DIR, exist_ok=True)

GTF = os.path.join(BASE, "stringtie_v4", "merged_stringtie_v4.gtf")
FA_NOVEL = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "novel_de_transcripts_only.fasta")

print("=" * 60)
print("  CTVT 36-sample — Step 2: Coding Potential Analysis")
print("=" * 60)

# 1. Load novel DE gene IDs
novel_de_genes = set()
with open(os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "novel_de_genes.txt")) as f:
    for line in f:
        novel_de_genes.add(line.strip())
print(f"  Novel DE genes: {len(novel_de_genes)}")

# 2. Build transcript → gene mapping from GTF
tx_to_gene = {}
gene_type = {}
with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != "transcript":
            continue
        attr = fields[8]
        gene_id = re.search(r'gene_id "([^"]+)"', attr)
        tx_id = re.search(r'transcript_id "([^"]+)"', attr)
        ref_gene = re.search(r'ref_gene_id "([^"]+)"', attr)
        gene_name = re.search(r'gene_name "([^"]+)"', attr)
        
        if gene_id and tx_id:
            gid = gene_id.group(1)
            tx_to_gene[tx_id.group(1)] = gid
            if gid not in gene_type:
                gene_type[gid] = "known" if (ref_gene or gene_name) else "novel"

print(f"  Transcripts in GTF: {len(tx_to_gene)}")

# 3. Identify novel DE transcript IDs
novel_de_tx_ids = set()
for tx_id, gene_id in tx_to_gene.items():
    if gene_id in novel_de_genes:
        novel_de_tx_ids.add(tx_id)

print(f"  Novel DE transcript IDs: {len(novel_de_tx_ids)}")

novel_seqs = {}
for record in SeqIO.parse(FA_NOVEL, "fasta"):
    tx_id = record.id.split()[0]
    if tx_id in novel_de_tx_ids:
        novel_seqs[tx_id] = record.seq

print(f"  Found {len(novel_seqs)} novel DE transcript sequences")

# 5. ORF analysis
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
                            aa = Seq(c).translate()
                            if str(aa) == '*':
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
    
    # Nucleotide composition
    n = len(seq)
    freqs = {'A': seq.count('A')/n, 'T': seq.count('T')/n, 
             'G': seq.count('G')/n, 'C': seq.count('C')/n}
    
    # Position asymmetry
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

# 6. Analyze all novel DE transcripts
print(f"\n  Analyzing coding potential of {len(novel_seqs)} transcripts...")
results = []

for tx_id, seq in novel_seqs.items():
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
    
    gene_id = tx_to_gene.get(tx_id, "")
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

df = pd.DataFrame(results)
df.to_csv(os.path.join(OUT_DIR, "coding_potential_results.tsv"), sep='\t', index=False)

# Split by classification
for cls in ['high_confidence_lncRNA', 'likely_lncRNA', 'ambiguous', 'potential_coding']:
    sub = df[df['classification'] == cls]
    sub.to_csv(os.path.join(OUT_DIR, f"{cls}.tsv"), sep='\t', index=False)
    print(f"  {cls}: {len(sub)} transcripts")

# Save gene IDs
lncrna_genes = set(df[df['classification'].isin(['high_confidence_lncRNA', 'likely_lncRNA'])]['gene_id'])
with open(os.path.join(OUT_DIR, "lncRNA_candidate_gene_ids.txt"), 'w') as f:
    for g in sorted(lncrna_genes):
        f.write(g + '\n')

print(f"\n  Combined lncRNA candidates (HC + Likely): {len(df[df['classification'].isin(['high_confidence_lncRNA', 'likely_lncRNA'])])} transcripts")
print(f"  Unique lncRNA gene IDs: {len(lncrna_genes)}")
print("\nStep 2 complete!")