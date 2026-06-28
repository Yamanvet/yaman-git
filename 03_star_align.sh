#!/bin/bash
# CTVT STAR Alignment — 4 samples at a time (16 threads each = 64 threads total)
# V3 settings (Q30 fastp trimmed reads)

STAR=/home/vet/miniconda3/envs/rnaseqY/bin/STAR
IDX=/home/vet/CTVT_raw_fastq_36samples/star_index_new
TRIM=/home/vet/CTVT_raw_fastq_36samples/fastp_test_v3
OUT=/home/vet/CTVT_raw_fastq_36samples/star_aligned_v3
TMP=/mnt/data8tb/star_tmp_v3

mkdir -p $OUT $TMP

# Get list of samples needing alignment
SAMPLES=()
for f in $TRIM/*_R1_trimmed.fastq.gz; do
    S=$(basename $f | sed 's/_R1_trimmed.fastq.gz//')
    [ -f "$OUT/${S}/${S}_Log.final.out" ] && continue
    SAMPLES+=("$S")
done

TOTAL=${#SAMPLES[@]}
echo "Samples to align: $TOTAL"
echo "Processing 4 at a time (16 threads each = 64 total)"
echo ""

BATCH=0
for ((i=0; i<TOTAL; i+=4)); do
    BATCH=$((BATCH+1))
    END=$((i+4))
    [ $END -gt $TOTAL ] && END=$TOTAL
    
    echo "=== BATCH $BATCH: samples $((i+1))-$END of $TOTAL ==="
    
    PIDS=()
    for ((j=i; j<END; j++)); do
        S=${SAMPLES[$j]}
        mkdir -p $OUT/$S
        
        systemd-run --user --scope -p MemoryMax=120G \
            $STAR --genomeDir $IDX \
            --readFilesIn $TRIM/${S}_R1_trimmed.fastq.gz $TRIM/${S}_R2_trimmed.fastq.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix $OUT/${S}/${S}_ \
            --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
            --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 \
            --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
            --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
            --alignSJoverhangMin 8 --quantMode GeneCounts \
            --outTmpDir $TMP/${S} --runThreadN 16 &
        PIDS+=($!)
        echo "  Started: $S"
    done
    
    echo "  Waiting for batch to complete..."
    for pid in ${PIDS[@]}; do
        wait $pid 2>/dev/null
    done
    
    COMPLETED=0
    for ((j=i; j<END; j++)); do
        S=${SAMPLES[$j]}
        if [ -f "$OUT/${S}/${S}_Log.final.out" ]; then
            COMPLETED=$((COMPLETED+1))
        fi
    done
    echo "  ✅ Batch $BATCH done: $COMPLETED/4 samples"
    echo ""
done

echo ""
echo "=== ALL DONE ==="
TOTAL_DONE=$(ls $OUT/*/*_Log.final.out 2>/dev/null | wc -l)
echo "Completed: $TOTAL_DONE/$TOTAL samples"
