#!/usr/bin/env python3
"""
CTVT: lncRNA cis-neighbor functional analysis (36 samples)
For each significant lncRNA DEG, find nearest protein-coding genes within 100kb
Then do enrichment on those cis-target protein-coding genes
"""

import csv, os, sys, collections

# ── 1. Load gene coordinates from merged GTF ──────────────────────────
print("Loading gene coordinates from GTF...")

import re
gene_coords = {}  # gene_id → (chr, start, end, strand, gene_name, gene_type)

gtf_file = "/home/vet/CTVT_raw_fastq_36samples/stringtie_v4/merged_stringtie_v4.gtf"
with open(gtf_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'transcript':
            continue
        
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        
        gene_id_match = re.search(r'gene_id "([^"]+)"', attrs)
        gene_name_match = re.search(r'gene_name "([^"]+)"', attrs)
        ref_gene_match = re.search(r'ref_gene_id "([^"]+)"', attrs)
        
        if not gene_id_match:
            continue
        
        gid = gene_id_match.group(1)
        gname = gene_name_match.group(1) if gene_name_match else (ref_gene_match.group(1) if ref_gene_match else gid)
        
        if gid not in gene_coords:
            gene_coords[gid] = (chrom, start, end, strand, gname)

print(f"  Loaded {len(gene_coords)} gene coordinates")

# ── 2. Load biotype classification ──────────────────────────────────
biotype_map = {}
for bt_file, bt_name in [
    ("/home/vet/CTVT_raw_fastq_36samples/tx2gene_protein_coding.csv", "protein_coding"),
    ("/home/vet/CTVT_raw_fastq_36samples/tx2gene_lncrna_known.csv", "lncRNA_known"),
]:
    with open(bt_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row['GENEID']
            if gid not in biotype_map:
                biotype_map[gid] = bt_name

# Mark remaining as novel
for gid in gene_coords:
    if gid not in biotype_map:
        biotype_map[gid] = 'novel'

# ── 3. Build sorted gene list per chromosome for fast neighbor search ──
from bisect import bisect_left, bisect_right

chrom_genes = collections.defaultdict(list)  # chr → [(start, end, gene_id)]
for gid, (chrom, start, end, strand, gname) in gene_coords.items():
    chrom_genes[chrom].append((start, end, gid))

for chrom in chrom_genes:
    chrom_genes[chrom].sort()

def find_neighbors(query_gid, window=100000):
    """Find protein-coding genes within window bp of query gene"""
    if query_gid not in gene_coords:
        return []
    
    qchrom, qstart, qend, qstrand, qgname = gene_coords[query_gid]
    
    neighbors = []
    if qchrom not in chrom_genes:
        return []
    
    genes = chrom_genes[qchrom]
    # Binary search for start position
    lo = bisect_left(genes, (qstart - window,))
    hi = bisect_right(genes, (qend + window + 1,))
    
    for gstart, gend, gid in genes[lo:hi+1]:
        if gid == query_gid:
            continue
        # Check distance
        if gstart > qend:
            dist = gstart - qend
        elif gend < qstart:
            dist = qstart - gend
        else:
            dist = 0  # overlapping
        
        if dist <= window and biotype_map.get(gid) == 'protein_coding':
            neighbors.append((gid, gene_coords[gid][4], dist, gene_coords[gid][0], gene_coords[gid][1], gene_coords[gid][2], gene_coords[gid][3]))
    
    neighbors.sort(key=lambda x: x[2])
    return neighbors

# ── 4. Load lncRNA DEGs and find cis-targets ──────────────────────────
base_dir = "/home/vet/CTVT_raw_fastq_36samples/deseq2_biotype_split_36/lncrna_known"
contrasts = ["Day2_vs_Day0", "Day6_vs_Day0", "Day6_vs_Day2", "Recovered_vs_Day0", "Recovered_vs_Day6"]

# Also load MSTRG → gene_name mapping
mstrg_to_name = {}
with open("/home/vet/CTVT_raw_fastq_36samples/mstrg_to_gene_name.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        mstrg_to_name[row['MSTRG_ID']] = row['gene_name']

all_results = {}

for cname in contrasts:
    print(f"\n{'='*60}")
    print(f"  {cname}")
    print(f"{'='*60}")
    
    deg_file = os.path.join(base_dir, f"{cname}_DEGs.csv")
    if not os.path.exists(deg_file):
        print("  File not found, skipping")
        continue
    
    degs = []
    with open(deg_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row['gene'].strip('"') if row.get('gene') else ''
            # Also try row name
            if not gid:
                gid = row.get('', '')
            degs.append(row)
    
    # Find cis targets for each lncRNA DEG
    cis_file = os.path.join(base_dir, f"{cname}_lncrna_cis_targets.csv")
    cis_rows = []
    all_cis_genes = set()
    
    for row in degs:
        gid = row.get('gene', row.get('', '')).strip('"')
        if not gid:
            # Try row name from first column
            gid = list(row.values())[0].strip('"')
        
        lncrna_name = mstrg_to_name.get(gid, gid)
        
        neighbors = find_neighbors(gid)
        for ngid, nname, dist, nchrom, nstart, nend, nstrand in neighbors:
            lfc = float(row.get('log2FoldChange', row.get('shrunkenLog2FC', '0')))
            padj = float(row.get('padj', '1'))
            cis_rows.append({
                'lncRNA_id': gid,
                'lncRNA_name': lncrna_name,
                'lncRNA_lfc': lfc,
                'lncRNA_padj': padj,
                'cis_gene_id': ngid,
                'cis_gene_name': nname,
                'distance': dist,
                'cis_chrom': nchrom,
                'cis_start': nstart,
                'cis_end': nend,
                'cis_strand': nstrand
            })
            all_cis_genes.add(nname)
    
    # Save cis targets
    with open(cis_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['lncRNA_id', 'lncRNA_name', 'lncRNA_lfc', 'lncRNA_padj',
                                                 'cis_gene_id', 'cis_gene_name', 'distance',
                                                 'cis_chrom', 'cis_start', 'cis_end', 'cis_strand'])
        writer.writeheader()
        writer.writerows(cis_rows)
    
    print(f"  {len(degs)} lncRNA DEGs → {len(cis_rows)} cis-target pairs → {len(all_cis_genes)} unique cis genes")

# ── 5. Summary ──────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("  SUMMARY")
print(f"{'='*60}")

for cname in contrasts:
    cis_file = os.path.join(base_dir, f"{cname}_lncrna_cis_targets.csv")
    deg_file = os.path.join(base_dir, f"{cname}_DEGs.csv")
    
    if os.path.exists(deg_file):
        with open(deg_file) as f:
            ndegs = sum(1 for _ in f) - 1
    else:
        ndegs = 0
    
    if os.path.exists(cis_file):
        with open(cis_file) as f:
            ncis = sum(1 for _ in f) - 1
        # Count unique genes
        genes = set()
        reader = csv.DictReader(open(cis_file))
        for row in reader:
            genes.add(row['cis_gene_name'])
        print(f"  {cname}: {ndegs} lncRNA DEGs → {ncis} cis pairs → {len(genes)} unique protein-coding targets")
    else:
        print(f"  {cname}: no cis targets file")

print("\nDone! cis-target files saved to lncrna_known/ directories")