#!/bin/bash
# ═══════════════════════════════════════════════════════════
# CTVT Phase 2: Salmon Quant — QC Test First
# 
# Runs salmon quant on 2 test samples to verify:
# - Mapping rate (should be 70-85% vs old 9.2%)
# - ISR strandedness flag correct
# - Gentrome index working
# Then proceeds to all 36 samples on confirmation
# ═══════════════════════════════════════════════════════════

set -euo pipefail

SSD_BASE="/home/vet/CTVT_raw_fastq_36samples"
FASTQ_DIR="$SSD_BASE/fastp_test_v3"
INDEX_DIR="$SSD_BASE/salmon_index/salmon_index_decoy"
QC_DIR="/mnt/data8tb/salmon_quant_QC"
FINAL_DIR="/mnt/data8tb/salmon_quant_final"
THREADS=16
LOG="$SSD_BASE/salmon_qc.log"

mkdir -p "$QC_DIR" "$FINAL_DIR"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

log "════════════════════════════════════════════════════"
log "  CTVT Phase 2: Salmon Quant — QC Test"
log "  Index: $INDEX_DIR"
log "  Fastq: $FASTQ_DIR"
log "════════════════════════════════════════════════════"

if [ ! -d "$INDEX_DIR" ]; then
    log "❌ Salmon index not found. Run gentrome_build.sh first."
    exit 1
fi

# ── QC Test: 2 samples ────────────────────────────────
log ""
log "══ QC TEST: 2 samples ══"
TEST_SAMPLES=("T1A_S1" "T1B_S2")

for S in "${TEST_SAMPLES[@]}"; do
    r1="$FASTQ_DIR/${S}_R1_trimmed.fastq.gz"
    r2="$FASTQ_DIR/${S}_R2_trimmed.fastq.gz"
    
    if [ ! -f "$r1" ]; then
        log "  ❌ FASTQ not found: $r1"
        continue
    fi
    
    log "  Quantifying $S..."
    salmon quant \
        -i "$INDEX_DIR" \
        -l ISR \
        -1 "$r1" -2 "$r2" \
        -o "$QC_DIR/${S}" \
        -p "$THREADS" \
        --validateMappings \
        --gcBias \
        >> "$LOG" 2>&1
    log "  ✅ $S done"
done

# ── Show mapping rates ────────────────────────────────
log ""
log "══ QC MAPPING RATES ══"
grep "Mapping rate" "$QC_DIR"/*/logs/salmon_quant.log 2>/dev/null || log "  (rates not found)"

log ""
log "══ QC COMPLETE ══"
log "Check mapping rates above."
log "If ≥70%, confirm to run all 36 samples."
log "Command: bash salmon_quant_full.sh"
