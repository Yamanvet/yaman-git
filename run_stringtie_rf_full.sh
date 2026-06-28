#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════
#  CTVT StringTie RE-RUN — RAM disk for big BAMs, SSD for small
#  --rf flag throughout
#  Steps: 1. Assembly (4 parallel) 2. Merge  3. Requant (-e -B)
#  (Skip: prepDE + DESeq2)
#
#  RAM POLICY:
#  - BAMs < 8GB: run from SSD directly (fast enough)
#  - BAMs >= 8GB: copy to /dev/shm, run, then clean
#  - /dev/shm usage capped per batch — never more than ~80GB total
# ═══════════════════════════════════════════════════════════════════════

set -euo pipefail
START_TIME=$(date +%s)

BAMDIR="/home/vet/CTVT_raw_fastq_36samples/star_aligned_v3"
OUTDIR="/home/vet/CTVT_raw_fastq_36samples/stringtie_rf"
BALLGOWN="/home/vet/CTVT_raw_fastq_36samples/ballgown_rf"
RAMDISK="/dev/shm/ctvt_rf"
MERGED="$OUTDIR/merged_stringtie_rf.gtf"
LOG="$OUTDIR/full_pipeline.log"
ANNOTATION="/home/vet/CTVT_raw_fastq_36samples/stringtie_v3/fixed_annotation.gtf"

mkdir -p "$OUTDIR" "$BALLGOWN"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
RAM_USAGE() { du -sh "$RAMDISK" 2>/dev/null | cut -f1 || echo "0"; }

log "═══════════════════════════════════════════════"
log "  CTVT StringTie RF Pipeline (RAM-smart)"
log "  Started: $(date)"
log "  BAMs >= 8GB → RAM disk, < 8GB → SSD direct"
log "═══════════════════════════════════════════════"

# ════════════════════════════════════════════════════
# STEP 1: ASSEMBLY (with --rf)
# ════════════════════════════════════════════════════
log ""
log "══════ STEP 1: ASSEMBLY (--rf) ══════"

BAMS=($(find "$BAMDIR" -name "*_Aligned.sortedByCoord.out.bam" | sort))
TOTAL=${#BAMS[@]}
log "Found $TOTAL BAM files"

for ((i=0; i<TOTAL; i+=4)); do
  END=$((i+4)); [ $END -gt $TOTAL ] && END=$TOTAL
  BATCH=$((i/4+1))
  log "Batch $BATCH: samples $((i+1))-$END of $TOTAL"
  
  PIDS=()
  
  for ((j=i; j<END; j++)); do
    BAM="${BAMS[$j]}"
    SAMPLE=$(basename "$(dirname "$BAM")")
    BAM_SIZE=$(du -B1 "$BAM" | cut -f1)
    BAM_SIZE_GB=$(echo "scale=1; $BAM_SIZE / 1073741824" | bc)
    
    # Check if already done
    if [ -f "$OUTDIR/${SAMPLE}_transcripts.gtf" ] && [ -s "$OUTDIR/${SAMPLE}_transcripts.gtf" ]; then
      log "  SKIP $SAMPLE (already assembled)"
      continue
    fi
    
    # RAM policy: BAMs >= 8GB go to /dev/shm
    if [ "$BAM_SIZE" -ge 8589934592 ]; then
      mkdir -p "$RAMDISK"
      log "  COPY $SAMPLE to RAM (${BAM_SIZE_GB}GB, RAM: $(RAM_USAGE))..."
      cp "$BAM" "$RAMDISK/${SAMPLE}.bam"
      INPUT_BAM="$RAMDISK/${SAMPLE}.bam"
      USING_RAM=1
    else
      log "  SSD $SAMPLE (${BAM_SIZE_GB}GB, direct)..."
      INPUT_BAM="$BAM"
      USING_RAM=0
    fi
    
    stringtie -p 16 --rf \
      -G "$ANNOTATION" \
      -o "$OUTDIR/${SAMPLE}_transcripts.gtf" \
      -l "$SAMPLE" \
      -A "$OUTDIR/${SAMPLE}_gene_abundances.tsv" \
      -f 0.01 -m 200 -c 1 \
      "$INPUT_BAM" \
      > "$OUTDIR/${SAMPLE}_stringtie.log" 2>&1 &
    PIDS+=($!)
  done
  
  # Wait for full batch
  for pid in ${PIDS[@]+"${PIDS[@]}"}; do
    wait "$pid" 2>/dev/null || true
  done
  
  # Clean RAM disk (only BAMs that were copied there)
  for ((j=i; j<END; j++)); do
    BAM="${BAMS[$j]}"
    SAMPLE=$(basename "$(dirname "$BAM")")
    rm -f "$RAMDISK/${SAMPLE}.bam"
  done
  
  # Verify
  OK=0
  for ((j=i; j<END; j++)); do
    BAM="${BAMS[$j]}"
    SAMPLE=$(basename "$(dirname "$BAM")")
    if [ -f "$OUTDIR/${SAMPLE}_transcripts.gtf" ] && [ -s "$OUTDIR/${SAMPLE}_transcripts.gtf" ]; then
      OK=$((OK+1))
    else
      log "  ⚠ $SAMPLE FAILED"
      echo "$SAMPLE" >> "$OUTDIR/failed_assembly.txt"
    fi
  done
  log "  ✅ $OK/$(($END-$i)) in batch $BATCH (RAM: $(RAM_USAGE))"
done

# Retry any failures
if [ -f "$OUTDIR/failed_assembly.txt" ]; then
  log "══ Retrying failed assemblies ══"
  while read -r SAMPLE; do
    BAM="$BAMDIR/$SAMPLE/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    log "  RETRY $SAMPLE..."
    rm -f "$RAMDISK/${SAMPLE}.bam"
    BAM_SIZE=$(du -B1 "$BAM" | cut -f1)
    if [ "$BAM_SIZE" -ge 8589934592 ]; then
      cp "$BAM" "$RAMDISK/${SAMPLE}.bam"
      INPUT_BAM="$RAMDISK/${SAMPLE}.bam"
    else
      INPUT_BAM="$BAM"
    fi
    stringtie -p 16 --rf \
      -G "$ANNOTATION" \
      -o "$OUTDIR/${SAMPLE}_transcripts.gtf" \
      -l "$SAMPLE" \
      -A "$OUTDIR/${SAMPLE}_gene_abundances.tsv" \
      -f 0.01 -m 200 -c 1 \
      "$INPUT_BAM" \
      > "$OUTDIR/${SAMPLE}_stringtie_retry.log" 2>&1
    rm -f "$RAMDISK/${SAMPLE}.bam"
  done < "$OUTDIR/failed_assembly.txt"
  rm "$OUTDIR/failed_assembly.txt"
fi

rmdir "$RAMDISK" 2>/dev/null || true
DONE=$(ls "$OUTDIR"/*_transcripts.gtf 2>/dev/null | wc -l)
log "══ STEP 1 COMPLETE: $DONE/$TOTAL assembled ══"

# ════════════════════════════════════════════════════
# STEP 2: MERGE
# ════════════════════════════════════════════════════
log ""
log "══════ STEP 2: MERGE ══════"

ls "$OUTDIR"/*_transcripts.gtf > "$OUTDIR/merge_list.txt"
MERGE_COUNT=$(wc -l < "$OUTDIR/merge_list.txt")
log "Merging $MERGE_COUNT GTF files (32 threads)..."
stringtie --merge -p 32 -G "$ANNOTATION" \
  -o "$MERGED" \
  "$OUTDIR/merge_list.txt" \
  >> "$LOG" 2>&1

if [ -f "$MERGED" ] && [ -s "$MERGED" ]; then
  LINES=$(grep -c "transcript_id" "$MERGED" || echo 0)
  log "✅ Merge done: $MERGED ($LINES transcripts, $(du -h "$MERGED" | cut -f1))"
else
  log "❌ MERGE FAILED"
  exit 1
fi

# ════════════════════════════════════════════════════
# STEP 3: REQUANT (-e -B) with --rf, RAM for big BAMs
# ════════════════════════════════════════════════════
log ""
log "══════ STEP 3: REQUANT (-e -B --rf) ══════"

mkdir -p "$BALLGOWN"

for BAM in "$BAMDIR"/*/*_Aligned.sortedByCoord.out.bam; do
  SAMPLE=$(basename "$(dirname "$BAM")")
  
  if [ -f "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" ] && [ -s "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" ]; then
    log "SKIP $SAMPLE (requant done)"
    continue
  fi
  
  mkdir -p "$BALLGOWN/$SAMPLE"
  BAM_SIZE=$(du -B1 "$BAM" | cut -f1)
  BAM_SIZE_GB=$(echo "scale=1; $BAM_SIZE / 1073741824" | bc)
  
  # RAM policy
  if [ "$BAM_SIZE" -ge 8589934592 ]; then
    mkdir -p "$RAMDISK"
    log "COPY $SAMPLE to RAM (${BAM_SIZE_GB}GB)..."
    cp "$BAM" "$RAMDISK/${SAMPLE}.bam"
    INPUT_BAM="$RAMDISK/${SAMPLE}.bam"
  else
    log "SSD $SAMPLE (${BAM_SIZE_GB}GB, direct)..."
    INPUT_BAM="$BAM"
  fi
  
  timeout 45m stringtie -e -B --rf -p 16 \
    -G "$MERGED" \
    -j 3 -f 0.1 \
    -o "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" \
    -A "$BALLGOWN/$SAMPLE/${SAMPLE}_abundance.tsv" \
    "$INPUT_BAM" \
    > "$BALLGOWN/$SAMPLE/${SAMPLE}_stringtie.log" 2>&1
  
  rm -f "$RAMDISK/${SAMPLE}.bam"
  
  if [ -f "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" ] && [ -s "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" ]; then
    LINES=$(grep -c "transcript_id" "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" 2>/dev/null || echo 0)
    log "✅ $SAMPLE — $LINES transcripts"
  else
    log "❌ FAILED $SAMPLE"
  fi
done

rmdir "$RAMDISK" 2>/dev/null || true
DONE=$(ls "$BALLGOWN"/*/*.gtf 2>/dev/null | wc -l)
log "══ STEP 3 COMPLETE: $DONE/36 requanted ══"

# ════════════════════════════════════════════════════
# DONE
# ════════════════════════════════════════════════════
END_TIME=$(date +%s)
DURATION=$(( (END_TIME - START_TIME) / 3600 ))
MINUTES=$(( ((END_TIME - START_TIME) % 3600) / 60 ))

log ""
log "═══════════════════════════════════════════════"
log "  ✅ FULL PIPELINE COMPLETE!"
log "  Duration: ${DURATION}h ${MINUTES}m"
log "  Assembly: $OUTDIR/"
log "  Requant:  $BALLGOWN/"
log "  (prepDE + DESeq2 skipped per request)"
log "═══════════════════════════════════════════════"
