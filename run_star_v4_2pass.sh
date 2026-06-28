#!/bin/bash
# ====================================================================
# CTVT: STAR v4 — 2-pass alignment with chimera detection + MD tags
#
# Runs 4 samples in parallel (4 × 16 threads = 64 total)
# Each process memory-limited: systemd-run --scope -p MemoryMax=120G
# Output + temp to 8TB drive
# ====================================================================

source /home/vet/miniconda3/etc/profile.d/conda.sh
conda activate rnaseqY

BASE="/home/vet/CTVT_raw_fastq_36samples"
FASTQ_DIR="$BASE/fastp_test_v3"
STAR_INDEX="$BASE/star_index_new"
OUTDIR="/mnt/data8tb/star_alignments"
TMPDIR="/mnt/data8tb/star_tmp_v3"
THREADS=16
MAX_PARALLEL=4

mkdir -p "$OUTDIR" "$TMPDIR"

# Get sample list
samples=$(ls "$FASTQ_DIR"/*_R1_trimmed.fastq.gz | sed 's/.*\///; s/_R1_trimmed.fastq.gz//' | sort)
total=$(echo "$samples" | wc -l)

echo "=========================================="
echo "  CTVT STAR v4 — 2-pass Alignment"
echo "  Samples: $total"
echo "  Parallel: $MAX_PARALLEL × ${THREADS} threads"
echo "  Output: $OUTDIR"
echo "=========================================="
echo ""

align_sample() {
    local sid=$1
    local count=$2
    local total=$3
    local r1="$FASTQ_DIR/${sid}_R1_trimmed.fastq.gz"
    local r2="$FASTQ_DIR/${sid}_R2_trimmed.fastq.gz"
    local outprefix="$OUTDIR/${sid}_"
    local logfile="$OUTDIR/${sid}_star.log"

    if [ -f "${outprefix}Aligned.sortedByCoord.out.bam" ]; then
        echo "[$count/$total] $sid — SKIP" >> "$OUTDIR/run_summary.log"
        return
    fi

    echo "[$count/$total] $sid — START $(date)" >> "$OUTDIR/run_summary.log"

    systemd-run --user --scope -p MemoryMax=120G -q \
    STAR --genomeDir "$STAR_INDEX" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$outprefix" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --outSAMstrandField intronMotif \
        --twopassMode Basic \
        --chimSegmentMin 12 \
        --chimOutType Junctions \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD \
        --quantMode GeneCounts \
        --outTmpDir "$TMPDIR/${sid}" \
        --runThreadN "$THREADS" \
        >> "$logfile" 2>&1

    echo "[$count/$total] $sid — DONE $(date)" >> "$OUTDIR/run_summary.log"
}

count=0
for sid in $samples; do
    count=$((count + 1))
    align_sample "$sid" "$count" "$total" &

    # Wait if we have MAX_PARALLEL jobs running
    running=$(jobs -r | wc -l)
    while [ "$running" -ge "$MAX_PARALLEL" ]; do
        sleep 10
        running=$(jobs -r | wc -l)
    done
done

# Wait for all remaining jobs
wait

echo ""
echo "=========================================="
echo "  Indexing BAMs..."
echo "=========================================="
for bam in "$OUTDIR"/*_Aligned.sortedByCoord.out.bam; do
    if [ -f "$bam" ] && [ ! -f "${bam}.bai" ]; then
        echo "Indexing: $(basename $bam)"
        samtools index "$bam"
    fi
done

echo ""
echo "=========================================="
echo "  Summary"
echo "=========================================="
echo "Total BAMs: $(ls "$OUTDIR"/*_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)/$total"
echo "Chimera junctions: $(ls "$OUTDIR"/*_Chimeric.out.junction 2>/dev/null | wc -l)"
echo "Done: $(date)"
echo "" >> "$OUTDIR/run_summary.log"
echo "All done: $(date)" >> "$OUTDIR/run_summary.log"
