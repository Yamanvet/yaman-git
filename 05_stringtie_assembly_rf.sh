#!/bin/bash
# CTVT StringTie Assembly with --rf (strand-specific reverse)
# ═══════════════════════════════════════════════════════════════════
# IMPORTANT: This version uses --rf flag for strand-specific RNA-seq
# The STAR alignment was run with --outSAMstrandField intronMotif
# which requires StringTie's --rf flag for proper transcript assembly
# ═══════════════════════════════════════════════════════════════════
#
# Results: +4.3% more transcripts vs no --rf (23.4M vs 22.4M total)
# 36/36 samples show improved transcript detection
#
OUTDIR="/home/vet/CTVT_raw_fastq_36samples/stringtie_rf"
FIXED="/home/vet/CTVT_raw_fastq_36samples/stringtie_v3/fixed_annotation.gtf"
BAMDIR="/home/vet/CTVT_raw_fastq_36samples/star_aligned_v3"

BAMS=($(find $BAMDIR -name "*_Aligned.sortedByCoord.out.bam" | sort))
TOTAL=${#BAMS[@]}
echo "Starting StringTie with --rf on $TOTAL BAMs (4 at a time, 16 threads each)..."

for ((i=0; i<TOTAL; i+=4)); do
  END=$((i+4)); [ $END -gt $TOTAL ] && END=$TOTAL
  BATCH=$((i/4+1))
  echo "[$(date +%H:%M)] Batch $BATCH: samples $((i+1))-$END of $TOTAL"
  
  PIDS=()
  for ((j=i; j<END; j++)); do
    SAMPLE=$(basename $(dirname "${BAMS[$j]}"))
    echo "  Starting $SAMPLE..."
    stringtie -p 16 --rf -G $FIXED \
      -o $OUTDIR/${SAMPLE}_transcripts.gtf \
      -l $SAMPLE \
      -A $OUTDIR/${SAMPLE}_gene_abundances.tsv \
      -f 0.01 -m 200 -c 1 \
      "${BAMS[$j]}" \
      > $OUTDIR/${SAMPLE}_stringtie.log 2>&1 &
    PIDS+=($!)
  done
  
  for pid in ${PIDS[@]}; do wait $pid 2>/dev/null; done
  
  OK=0
  for ((j=i; j<END; j++)); do
    SAMPLE=$(basename $(dirname "${BAMS[$j]}"))
    [ -s "$OUTDIR/${SAMPLE}_transcripts.gtf" ] && OK=$((OK+1)) || echo "  ⚠️ $SAMPLE FAILED"
  done
  echo "  ✅ $OK/$(($END-$i)) done"
done

DONE=$(ls $OUTDIR/*_transcripts.gtf 2>/dev/null | wc -l)
echo ""
echo "[$(date +%H:%M)] StringTie Assembly Complete: $DONE/$TOTAL"
