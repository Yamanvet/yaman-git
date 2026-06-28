#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════
#  CTVT StringTie v4 — New settings (RAM-smart + BAM index)
#  --rf throughout
#  Steps: 1. Assembly (4 parallel) 2. Merge (f=0.05) 3. Requant (-e -B)
#
#  Key changes from v3:
#  - New params: -j 5 -c 2.5 -f 0.1 -m 200  (vs old -f 0.01 -m 200 -c 1)
#  - Merge param: -f 0.05 (was default)
#  - BAMs have .bai index (faster StringTie lookup)
#  - RAM: BAMs >= 8GB → /dev/shm, rest → direct
# ═══════════════════════════════════════════════════════════════════════

set -euo pipefail
START_TIME=$(date +%s)

BAMDIR="/mnt/data8tb/star_alignments"
OUTDIR="/home/vet/CTVT_raw_fastq_36samples/stringtie_v4"
BALLGOWN="/home/vet/CTVT_raw_fastq_36samples/ballgown_v4"
RAMDISK="/dev/shm/ctvt_v4"
MERGED="$OUTDIR/merged_stringtie_v4.gtf"
LOG="$OUTDIR/full_pipeline.log"
ANNOTATION="/mnt/data8tb/CTVT_Raw_Data_root_backup/clean_stringtie_annotation.gtf"

mkdir -p "$OUTDIR" "$BALLGOWN"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
RAM_USAGE() { du -sh "$RAMDISK" 2>/dev/null | cut -f1 || echo "0"; }

log "═══════════════════════════════════════════════"
log "  CTVT StringTie v4 Pipeline"
log "  Started: $(date)"
log "  Settings: -j 5 -c 2.5 -f 0.1 (assembly), -f 0.05 (merge)"
log "  BAMs >= 8GB → RAM disk, < 8GB → SSD direct"
log "═══════════════════════════════════════════════"

# ════════════════════════════════════════════════════
# STEP 1: ASSEMBLY (with --rf, new params)
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
    SAMPLE=$(basename "$BAM" | sed 's/_Aligned.sortedByCoord.out.bam//')
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
      cp "${BAM}.bai" "$RAMDISK/${SAMPLE}.bam.bai" 2>/dev/null || true
      INPUT_BAM="$RAMDISK/${SAMPLE}.bam"
    else
      log "  SSD $SAMPLE (${BAM_SIZE_GB}GB, direct)..."
      INPUT_BAM="$BAM"
    fi
    
    stringtie -p 8 --rf \
      -G "$ANNOTATION" \
      -o "$OUTDIR/${SAMPLE}_transcripts.gtf" \
      -l "$SAMPLE" \
      -A "$OUTDIR/${SAMPLE}_gene_abundances.tsv" \
      -j 5 -c 2.5 -f 0.1 -m 200 \
      "$INPUT_BAM" \
      > "$OUTDIR/${SAMPLE}_stringtie.log" 2>&1 &
    PIDS+=($!)
  done
  
  # Wait for full batch
  for pid in ${PIDS[@]+"${PIDS[@]}"}; do
    wait "$pid" 2>/dev/null || true
  done
  
  # Clean RAM disk
  for ((j=i; j<END; j++)); do
    BAM="${BAMS[$j]}"
    SAMPLE=$(basename "$BAM" | sed 's/_Aligned.sortedByCoord.out.bam//')
    rm -f "$RAMDISK/${SAMPLE}.bam" "$RAMDISK/${SAMPLE}.bam.bai"
  done
  
  # Verify
  OK=0
  for ((j=i; j<END; j++)); do
    BAM="${BAMS[$j]}"
    SAMPLE=$(basename "$BAM" | sed 's/_Aligned.sortedByCoord.out.bam//')
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
    BAM="$BAMDIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    log "  RETRY $SAMPLE..."
    rm -f "$RAMDISK/${SAMPLE}.bam" "$RAMDISK/${SAMPLE}.bam.bai"
    BAM_SIZE=$(du -B1 "$BAM" | cut -f1)
    if [ "$BAM_SIZE" -ge 8589934592 ]; then
      cp "$BAM" "$RAMDISK/${SAMPLE}.bam"
      cp "${BAM}.bai" "$RAMDISK/${SAMPLE}.bam.bai" 2>/dev/null || true
      INPUT_BAM="$RAMDISK/${SAMPLE}.bam"
    else
      INPUT_BAM="$BAM"
    fi
    stringtie -p 8 --rf \
      -G "$ANNOTATION" \
      -o "$OUTDIR/${SAMPLE}_transcripts.gtf" \
      -l "$SAMPLE" \
      -A "$OUTDIR/${SAMPLE}_gene_abundances.tsv" \
      -j 5 -c 2.5 -f 0.1 -m 200 \
      "$INPUT_BAM" \
      > "$OUTDIR/${SAMPLE}_stringtie_retry.log" 2>&1
    rm -f "$RAMDISK/${SAMPLE}.bam" "$RAMDISK/${SAMPLE}.bam.bai"
  done < "$OUTDIR/failed_assembly.txt"
  rm "$OUTDIR/failed_assembly.txt"
fi

rmdir "$RAMDISK" 2>/dev/null || true
DONE=$(ls "$OUTDIR"/*_transcripts.gtf 2>/dev/null | wc -l)
log "══ STEP 1 COMPLETE: $DONE/$TOTAL assembled ══"

# ════════════════════════════════════════════════════
# STEP 2: MERGE (with -f 0.05)
# ════════════════════════════════════════════════════
log ""
log "══════ STEP 2: MERGE ══════"

ls "$OUTDIR"/*_transcripts.gtf > "$OUTDIR/merge_list.txt"
MERGE_COUNT=$(wc -l < "$OUTDIR/merge_list.txt")
log "Merging $MERGE_COUNT GTF files (16 threads)..."
stringtie --merge -p 16 \
  -G "$ANNOTATION" \
  -f 0.05 -m 200 \
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
# STEP 3: REQUANT (-e -B --rf) with RAM optimization
# ════════════════════════════════════════════════════
log ""
log "══════ STEP 3: REQUANT (-e -B --rf) ══════"

for BAM in "$BAMDIR"/*_Aligned.sortedByCoord.out.bam; do
  SAMPLE=$(basename "$BAM" | sed 's/_Aligned.sortedByCoord.out.bam//')
  
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
    cp "${BAM}.bai" "$RAMDISK/${SAMPLE}.bam.bai" 2>/dev/null || true
    INPUT_BAM="$RAMDISK/${SAMPLE}.bam"
  else
    log "SSD $SAMPLE (${BAM_SIZE_GB}GB, direct)..."
    INPUT_BAM="$BAM"
  fi
  
  timeout 45m stringtie -e -B --rf -p 8 \
    -G "$MERGED" \
    -o "$BALLGOWN/$SAMPLE/${SAMPLE}.gtf" \
    -A "$BALLGOWN/$SAMPLE/${SAMPLE}_abundance.tsv" \
    "$INPUT_BAM" \
    > "$BALLGOWN/$SAMPLE/${SAMPLE}_stringtie.log" 2>&1
  
  rm -f "$RAMDISK/${SAMPLE}.bam" "$RAMDISK/${SAMPLE}.bam.bai"
  
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
log "═══════════════════════════════════════════════"
