#!/bin/bash
set -euo pipefail

###########################################################
# IQ-TREE: Maximum-likelihood phylogeny with validation
# Input: Gubbins filtered SNP alignment
###########################################################

IN="gubbins/campy.filtered_polymorphic_sites.fasta"
OUT_DIR="iqtree"
LOG_DIR="logs"
PREFIX="campy"
THREADS=63   # <— set your number of CPU threads here

mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "==========================================="
echo "[IQ-TREE] Starting ML phylogeny inference"
echo "Input alignment: $IN"
echo "Output prefix:   $OUT_DIR/$PREFIX"
echo "Threads:         $THREADS"
echo "==========================================="

# 1. Check input exists
if [[ ! -f "$IN" ]]; then
    echo "[ERROR] Input alignment not found: $IN"
    exit 1
fi

# 2. Run IQ-TREE
iqtree2 \
    -s "$IN" \
    -m GTR+F+I+G4 \
    -bb 1000 \
    -nt "$THREADS" \
    -pre "$OUT_DIR/$PREFIX" \
    > "$LOG_DIR/iqtree.log" 2>&1

# 3. Check for errors in log
if grep -qi "error" "$LOG_DIR/iqtree.log"; then
    echo "[ERROR] IQ-TREE reported an error. Check logs/iqtree.log"
    exit 1
fi

# 4. Check expected output exists
TREE="${OUT_DIR}/${PREFIX}.treefile"
IQLOG="${OUT_DIR}/${PREFIX}.iqtree"

if [[ ! -f "$TREE" ]]; then
    echo "[ERROR] Tree file missing: $TREE"
    exit 1
fi

if [[ ! -f "$IQLOG" ]]; then
    echo "[ERROR] IQ-TREE statistics file missing: $IQLOG"
    exit 1
fi

# 5. Check tree is non-empty
if [[ ! -s "$TREE" ]]; then
    echo "[ERROR] Tree file is empty: $TREE"
    exit 1
fi

echo "==========================================="
echo "[IQ-TREE] Success!"
echo "Final ML tree: $TREE"
echo "IQ-TREE report: $IQLOG"
echo "==========================================="
