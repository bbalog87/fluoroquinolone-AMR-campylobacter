#!/bin/bash
set -euo pipefail

###############################################
# Gubbins: recombination filtering with validation
###############################################

IN="panaroo/core_gene_alignment_filtered.aln"
OUT_DIR="gubbins"
LOG_DIR="logs"
PREFIX="campy"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# 1. Check input exists
if [[ ! -f "$IN" ]]; then
    echo "[ERROR] Input alignment not found: $IN"
    exit 1
fi

echo "[GUBBINS] Starting recombination filtering"
echo "Input:  $IN"
echo "Output: $OUT_DIR/$PREFIX.*"

# 2. Run Gubbins
run_gubbins.py \
    --prefix "$OUT_DIR/$PREFIX" \
    --threads 63 \
    "$IN" \
    > "$LOG_DIR/gubbins.log" 2>&1

# 3. Check for errors in log
if grep -qi "error" "$LOG_DIR/gubbins.log"; then
    echo "[ERROR] Gubbins reported an error. Check logs/gubbins.log"
    exit 1
fi

# 4. Check expected output exists
FILTERED="${OUT_DIR}/${PREFIX}.filtered_polymorphic_sites.fasta"
TREE="${OUT_DIR}/${PREFIX}.final_tree.tre"

if [[ ! -f "$FILTERED" ]]; then
    echo "[ERROR] Filtered SNP alignment missing: $FILTERED"
    exit 1
fi

if [[ ! -f "$TREE" ]]; then
    echo "[ERROR] Final tree missing: $TREE"
    exit 1
fi

# 5. Check filtered alignment is non-empty
if [[ ! -s "$FILTERED" ]]; then
    echo "[ERROR] Filtered SNP alignment is empty."
    exit 1
fi

echo "[GUBBINS] Success!"
echo "Filtered SNP alignment: $FILTERED"
echo "Final tree:             $TREE"
