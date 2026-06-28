#!/bin/bash
# ═══════════════════════════════════════════════════════════
# CTVT Phase 2b: Salmon Quant — All 36 Samples
# Run only after QC test confirms ≥70% mapping rate
# ═══════════════════════════════════════════════════════════

set -euo pipefail
START_TIME=$(date +%s)

SSD_BASE="/home/vet/CTVT_raw_fastq_36samples"
FASTQ_DIR="$SSD_BASE/fastp_test_v3"
INDEX_DIR="$SSD_BASE/salmon_index/salmon_index_decoy"
OUTDIR="/mnt/data8tb/salmon_quant_final"
LOG="$SSD_BASE/salmon_quant_full.log"
THREADS=16

mkdir -p "$OUTDIR"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

log "════════════════════════════════════════════════════"
log "  CTVT: Salmon Quant — All 36 Samples"
log "  Started: $(date)"
log "  Threads: $THREADS"
log "════════════════════════════════════════════════════"

SAMPLES=($(ls "$FASTQ_DIR"/*_R1_trimmed.fastq.gz | \
    sed 's|.*/||; s|_R1_trimmed.fastq.gz||' | sort))
TOTAL=${#SAMPLES[@]}
log "Found $TOTAL samples"

# Run 2 in parallel to balance speed vs resource usage
for ((i=0; i<TOTAL; i+=2)); do
    END=$((i+2)); [ $END -gt $TOTAL ] && END=$TOTAL
    BATCH=$((i/2+1))
    log "Batch $BATCH: samples $((i+1))-$END"
    
    for ((j=i; j<END; j++)); do
        S="${SAMPLES[$j]}"
        r1="$FASTQ_DIR/${S}_R1_trimmed.fastq.gz"
        r2="$FASTQ_DIR/${S}_R2_trimmed.fastq.gz"
        
        if [ -d "$OUTDIR/$S" ] && [ -f "$OUTDIR/$S/quant.sf" ]; then
            log "  SKIP $S (already done)"
            continue
        fi
        
        salmon quant \
            -i "$INDEX_DIR" \
            -l ISR \
            -1 "$r1" -2 "$r2" \
            -o "$OUTDIR/$S" \
            -p "$THREADS" \
            --validateMappings \
            --gcBias \
            >> "$LOG" 2>&1 &
    done
    
    wait
    log "  ✅ Batch $BATCH done"
done

# Summary
END_TIME=$(date +%s)
DURATION=$(( (END_TIME - START_TIME) / 60 ))
DONE=$(ls -d "$OUTDIR"/*/ 2>/dev/null | wc -l)

log ""
log "════════════════════════════════════════════════════"
log "  ✅ PHASE 2 COMPLETE"
log "  Duration: ${DURATION} min"
log "  Samples: $DONE/$TOTAL"
log "  Output: $OUTDIR"
log "════════════════════════════════════════════════════"
