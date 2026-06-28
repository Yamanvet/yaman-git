#!/bin/bash
# CTVT StringTie Re-quant — RAW BAM on RAM DISK (--rf strand-specific)
OUTDIR="/home/vet/CTVT_raw_fastq_36samples/ballgown_rf"
MERGED="/home/vet/CTVT_raw_fastq_36samples/stringtie_rf/merged_stringtie_v3.gtf"
BAMDIR="/home/vet/CTVT_raw_fastq_36samples/star_aligned_v3"
RAMDISK="/dev/shm/ctvt"
LOGFILE="$OUTDIR/requant_ram.log"

echo "==================================================================" | tee "$LOGFILE"
echo "$(date): PURE RAM re-quant (--rf strand-specific) — sequential" | tee -a "$LOGFILE"
echo "==================================================================" | tee -a "$LOGFILE"

for bam in "$BAMDIR"/*/*_Aligned.sortedByCoord.out.bam; do
  sample=$(basename "$(dirname "$bam")")
  
  if [ -f "$OUTDIR/$sample/${sample}.gtf" ] && [ -s "$OUTDIR/$sample/${sample}.gtf" ]; then
    echo "$(date): SKIP $sample (already done)" | tee -a "$LOGFILE"
    continue
  fi
  
  mkdir -p "$OUTDIR/$sample" "$RAMDISK"
  bam_mb=$(du -m "$bam" | cut -f1)
  
  # Copy to RAM
  echo "$(date): COPY $sample → /dev/shm ($bam_mb MB)..." | tee -a "$LOGFILE"
  cp "$bam" "$RAMDISK/${sample}.bam"
  echo "$(date): DONE copy $sample" | tee -a "$LOGFILE"
  
  # Run StringTie from RAM
  echo "$(date): STRINGTIE $sample..." | tee -a "$LOGFILE"
  stringtie -e -B --rf -p 16 -G "$MERGED" \
    -j 3 -f 0.1 \
    -o "$OUTDIR/${sample}/${sample}.gtf" \
    -A "$OUTDIR/${sample}/${sample}_abundance.tsv" \
    "$RAMDISK/${sample}.bam" > "$OUTDIR/${sample}/${sample}_stringtie.log" 2>&1
  
  # Clean RAM
  rm -f "$RAMDISK/${sample}.bam"
  
  if [ -f "$OUTDIR/$sample/${sample}.gtf" ] && [ -s "$OUTDIR/$sample/${sample}.gtf" ]; then
    n=$(grep -c "transcript_id" "$OUTDIR/$sample/${sample}.gtf" 2>/dev/null || echo 0)
    echo "$(date): ✅ $sample — $n transcripts" | tee -a "$LOGFILE"
  else
    echo "$(date): ❌ FAILED $sample" | tee -a "$LOGFILE"
  fi
  echo "" | tee -a "$LOGFILE"
done

rmdir "$RAMDISK" 2>/dev/null
completed=$(ls "$OUTDIR"/*/*.gtf 2>/dev/null | wc -l)
echo "$(date): RUN COMPLETE — $completed/36 done" | tee -a "$LOGFILE"
