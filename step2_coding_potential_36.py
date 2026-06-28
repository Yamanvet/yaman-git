#!/usr/bin/env python3
"""
CTVT 36-sample: Step 2 — Coding potential analysis for novel lncRNA candidates
ORF analysis + Fickett score + classification
"""
import os, re, csv
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from collections import defaultdict

BASE = "/home/vet/CTVT_raw_fastq_36samples"
OUT_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "coding_potential")
os.makedirs(OUT_DIR, exist_ok=True)

GTF = os.path.join(BASE, "stringtie_v4", "merged_stringtie_v4.gtf")
FA = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline", "DE_novel_transcripts.fasta")

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

# 3. ORF finder
def find_longest_orf(seq, min_length=30):
    """Find longest ORF in all 6 frames."""
    seq = str(seq).upper().replace('N', '')
    if len(seq) < 3:
        return None, 0, 0, 0
    
    best_orf = None
    best_len = 0
    best_frame = 0
    
    for strand, sequence in [(1, seq), (-1, str(Seq(seq).reverse_complement()))]:
        for frame in range(3):
            i = frame
            while i < len(sequence) - 2:
                codon = sequence[i:i+3]
                if codon == 'ATG':
                    # Found start codon, find ORF end
                    aa_len = 0
                    j = i
                    while j < len(sequence) - 2:
                        codon = sequence[j:j+3]
                        aa = Seq(codon).translate()
                        if aa == '*':
                            break
                        aa_len += 1
                        j += 3
                    if aa_len > best_len:
                        best_len = aa_len
                        best_frame = strand * (frame + 1)
                        best_orf = sequence[i:j+3]
                i += 3
    
    return best_orf, best_len, best_frame, len(best_orf) if best_orf else 0

# 4. Fickett score (simplified — nucleotide composition bias)
def fickett_score(seq):
    """Calculate Fickett score based on nucleotide composition."""
    seq = str(seq).upper()
    if len(seq) < 2:
        return 0.0
    
    A = seq.count('A') / len(seq)
    T = seq.count('T') / len(seq)
    G = seq.count('G') / len(seq)
    C = seq.count('C') / len(seq)
    
    # Position asymmetry
    posA = posT = posG = posC = 0
    for i, nt in enumerate(seq):
        if i % 3 == 0:
            if nt == 'A': posA += 1
            elif nt == 'T': posT += 1
            elif nt == 'G': posG += 1
            elif nt == 'C': posC += 1
    
    n_thirds = len(seq) // 3
    if n_thirds == 0:
        return 0.0
    
    posA /= n_thirds
    posT /= n_thirds
    posG /= n_thirds
    posC /= n_thirds
    
    # Simplified Fickett: higher for coding, lower for non-coding
    score = (A * posA + T * posT + G * posG + C * posC) / 4
    return min(score, 1.0)

# 5. Load FASTA and analyze
print("\n  Analyzing coding potential...")
results = []
seen = 0
novel_only = 0

for record in SeqIO.parse(FA, "fasta"):
    tx_id = record.id.split()[0]
    
    # Only analyze novel DE gene transcripts
    gene_id = tx_to_gene.get(tx_id, "")
    if gene_id not in novel_de_genes:
        continue
    
    seq = str(record.seq)
    seq_len = len(seq)
    
    # Skip very short sequences
    if seq_len < 30:
        continue
    
    # Find longest ORF
    orf, orf_aa_len, orf_frame, orf_nt_len = find_longest_orf(seq)
    
    # Calculate ORF coverage
    orf_coverage = orf_nt_len / seq_len if seq_len > 0 else 0
    
    # Fickett score
    fickett = fickett_score(seq)
    
    # Classification
    if orf_aa_len == 0 or orf_aa_len < 30:
        classification = "high_confidence_lncRNA"
    elif orf_aa_len < 100 and orf_coverage < 0.1 and fickett < 0.1:
        classification = "high_confidence_lncRNA"
    elif orf_aa_len < 100:
        classification = "likely_lncRNA"
    elif orf_aa_len < 300 and orf_coverage < 0.3:
        classification = "likely_lncRNA"
    elif orf_aa_len < 300:
        classification = "ambiguous"
    else:
        classification = "potential_coding"
    
    results.append({
        'transcript_id': tx_id,
        'gene_id': gene_id,
        'seq_length': seq_len,
        'longest_ORF_aa': orf_aa_len,
        'ORF_coverage': round(orf_coverage, 3),
        'ORF_frame': orf_frame,
        'fickett_score': round(fickett, 4),
        'classification': classification
    })
    
    novel_only += 1
    seen += 1
    
    if seen % 5000 == 0:
        print(f"    Processed {seen} novel DE transcripts...")

print(f"\n  Total novel DE transcripts analyzed: {len(results)}")

# 6. Save results
df = pd.DataFrame(results)  # noqa: imported at top
# Actually need to import pandas here
import pandas as pd
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

print(f"\n  Combined lncRNA candidates: {len(df[df['classification'].isin(['high_confidence_lncRNA', 'likely_lncRNA'])])} transcripts")
print(f"  Unique lncRNA gene IDs: {len(lncrna_genes)}")
print("\nStep 2 complete!")