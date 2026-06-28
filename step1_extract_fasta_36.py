#!/usr/bin/env python3
"""
CTVT 36-sample: Step 1 — Extract FASTA for novel DE transcripts
Uses gffread to extract from genome FASTA using the merged GTF.
"""
import os, re, subprocess
import pandas as pd

BASE = "/home/vet/CTVT_raw_fastq_36samples"
DEG_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "protein_coding")
GTF = os.path.join(BASE, "stringtie_v4", "merged_stringtie_v4.gtf")
GENOME_FA = "/mnt/data8tb/CTVT_Raw_Data_root_backup/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna"
OUT_DIR = os.path.join(BASE, "deseq2_biotype_split_36", "lncrna_pipeline")
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("  CTVT 36-sample — Step 1: FASTA Extraction")
print("=" * 60)

# 1. Collect all novel DE gene IDs from the novel lncRNA candidates
novel_genes = set()
contrasts = ["Day2_vs_Day0", "Day6_vs_Day0", "Day6_vs_Day2", "Recovered_vs_Day0", "Recovered_vs_Day6"]

for cname in contrasts:
    f = os.path.join(BASE, "deseq2_biotype_split_36", "novel_lncRNA_candidates", f"{cname}_novel_MSTRG.csv")
    if os.path.exists(f):
        df = pd.read_csv(f, index_col=0)
        novel_genes.update(df.index)

print(f"\n  Total unique novel MSTRG genes in DE results: {len(novel_genes)}")

# 2. Build transcript annotation from GTF
print("  Building gene annotation map from GTF...")
mstrg_to_name = {}
mstrg_to_ref = {}
gene_type_map = {}
transcript_to_gene = {}
transcript_info = {}  # tx_id -> (chr, start, end, strand, exon_count, tx_len)

with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        feature = fields[2]
        attr = fields[8]
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        
        gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
        if not gene_id_match:
            continue
        gene_id = gene_id_match.group(1)
        
        gene_name_match = re.search(r'gene_name "([^"]+)"', attr)
        ref_gene_match = re.search(r'ref_gene_id "([^"]+)"', attr)
        tx_id_match = re.search(r'transcript_id "([^"]+)"', attr)
        
        gene_name = gene_name_match.group(1) if gene_name_match else None
        ref_gene = ref_gene_match.group(1) if ref_gene_match else None
        
        if ref_gene or gene_name:
            mstrg_to_name[gene_id] = gene_name or ref_gene
            mstrg_to_ref[gene_id] = ref_gene
            gene_type_map[gene_id] = "known"
        else:
            gene_type_map[gene_id] = "novel"
        
        if feature == "transcript" and tx_id_match:
            tx_id = tx_id_match.group(1)
            transcript_to_gene[tx_id] = gene_id

# 3. Collect transcript IDs for novel DE genes
novel_tx_ids = set()
with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        if fields[2] != "transcript":
            continue
        attr = fields[8]
        gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
        tx_id_match = re.search(r'transcript_id "([^"]+)"', attr)
        if gene_id_match and tx_id_match:
            gid = gene_id_match.group(1)
            if gid in novel_genes:
                novel_tx_ids.add(tx_id_match.group(1))

print(f"  Novel DE gene IDs: {len(novel_genes)}")
print(f"  Corresponding transcript IDs: {len(novel_tx_ids)}")

# 4. Extract transcript sequences using gffread
print("\n  Extracting transcript sequences with gffread...")
out_fa = os.path.join(OUT_DIR, "DE_novel_transcripts.fasta")

cmd = f"gffread -g {GENOME_FA} -w {out_fa} {GTF}"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
if result.returncode != 0:
    print(f"  gffread error: {result.stderr[:200]}")
else:
    print(f"  Extracted transcripts to {out_fa}")

# 5. Create annotation summary
print("\n  Creating transcript annotation summary...")
annot_data = []
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
        gene_name = re.search(r'gene_name "([^"]+)"', attr)
        ref_gene = re.search(r'ref_gene_id "([^"]+)"', attr)
        
        if gene_id and tx_id:
            gid = gene_id.group(1)
            annot_data.append({
                'transcript_id': tx_id.group(1),
                'gene_id': gid,
                'gene_name': gene_name.group(1) if gene_name else (ref_gene.group(1) if ref_gene else gid),
                'ref_gene_id': ref_gene.group(1) if ref_gene else '',
                'gene_type': gene_type_map.get(gid, 'novel'),
                'is_DE_novel': gid in novel_genes
            })

annot_df = pd.DataFrame(annot_data)
annot_df.to_csv(os.path.join(OUT_DIR, "transcript_annotation_summary.tsv"), sep='\t')
print(f"  Saved transcript annotation: {len(annot_df)} transcripts")
print(f"    Known: {sum(annot_df['gene_type'] == 'known')}")
print(f"    Novel: {sum(annot_df['gene_type'] == 'novel')}")
print(f"    DE novel: {sum(annot_df['is_DE_novel'])}")

print("\nStep 1 complete!")