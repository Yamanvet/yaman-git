#!/bin/bash
# ====================================================================
# Salmon quantification — all 36 CTVT samples
# Uses trimmed FASTQs from fastp_test_v3
# ====================================================================

source /home/vet/miniconda3/etc/profile.d/conda.sh
conda activate deseq2_env

BASE="/home/vet/CTVT_raw_fastq_36samples"
FASTQ_DIR="$BASE/fastp_test_v3"
SALMON_INDEX="$BASE/salmon_index/canine_index"
QUANT_DIR="$BASE/salmon_quant"
mkdir -p "$QUANT_DIR"

# Get all sample IDs from R1 files
samples=$(ls "$FASTQ_DIR"/*_R1_trimmed.fastq.gz | sed 's/.*\///; s/_R1_trimmed.fastq.gz//')

echo "Starting Salmon quantification for $(echo "$samples" | wc -l) samples"
echo "Using ${SALMON_INDEX}"

count=0
total=$(echo "$samples" | wc -l)

for sid in $samples; do
    count=$((count + 1))
    r1="$FASTQ_DIR/${sid}_R1_trimmed.fastq.gz"
    r2="$FASTQ_DIR/${sid}_R2_trimmed.fastq.gz"
    out="$QUANT_DIR/$sid"
    
    if [ -f "$out/quant.sf" ]; then
        echo "[$count/$total] $sid — SKIP (already done)"
        continue
    fi
    
    echo "[$count/$total] $sid — quantifying..."
    salmon quant \
        -i "$SALMON_INDEX" \
        -l A \
        -1 "$r1" \
        -2 "$r2" \
        -o "$out" \
        --threads 4 \
        --validateMappings \
        --gcBias \
        2>&1 | grep -E "processed|mapping rate|num_compatible"
    
    echo "  Done: $sid"
done

echo ""
echo "All samples quantified!"
