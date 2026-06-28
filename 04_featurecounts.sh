#!/usr/bin/env bash
set -euo pipefail

# ========== V3 featureCounts ==========
GTF="/mnt/data8tb/CTVT_Raw_Data_root_backup/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gtf"
BAM_DIR="/home/vet/CTVT_raw_fastq_36samples/star_aligned_v3"
OUT_DIR="/home/vet/CTVT_raw_fastq_36samples/featureCounts_v3"
THREADS=16
# ======================================

mkdir -p "$OUT_DIR"

# collect BAMs
echo "Collecting BAMs from $BAM_DIR..."
BAM_LIST=$(find "$BAM_DIR" -name "*_Aligned.sortedByCoord.out.bam" | sort)

BAM_COUNT=$(echo "$BAM_LIST" | wc -l)
echo "Found $BAM_COUNT BAM files."

if [ "$BAM_COUNT" -eq 0 ]; then
  echo "ERROR: No BAMs found!" >&2
  exit 1
fi

echo "Starting featureCounts with -s 0 (unstranded)..."
echo ""
featureCounts -T "$THREADS" \
  -p --countReadPairs \
  -B -C \
  -t exon \
  -g gene_id \
  -s 0 \
  -a "$GTF" \
  -o "$OUT_DIR/gene_counts_v3.txt" \
  $BAM_LIST

echo ""
echo "Done! Output: $OUT_DIR/gene_counts_v3.txt"
echo "Summary:"
head -15 "$OUT_DIR/gene_counts_v3.txt.summary"
