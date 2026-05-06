#!/bin/bash
set -euo pipefail

###############################################
# Panaroo pangenome + core genome alignment
# Campylobacter phylogenomics pipeline
###############################################

INPUT_DIR="prokka"
OUTPUT_DIR="panaroo"
LOG_DIR="logs"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "[PANAROO] Running pangenome + core alignment"

panaroo \
    -i $INPUT_DIR/*/*.gff \
    -o $OUTPUT_DIR \
    --clean-mode strict \
    --remove-invalid-genes \
    --aligner mafft \
    --core_threshold 0.95 \
    --alignment core \
    --threads 64 \
    > $LOG_DIR/panaroo.log 2>&1

echo "[PANAROO] Finished"
echo "[PANAROO] Core alignment: $OUTPUT_DIR/core_alignment.aln"

