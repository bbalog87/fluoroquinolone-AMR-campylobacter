#!/bin/bash
set -euo pipefail

INPUT=genomes
OUT=prokka
LOG=logs

mkdir -p $OUT $LOG

for f in $INPUT/*.fna; do
    base=$(basename "$f" .fna)
    echo "[PROKKA] Annotating $base"
    prokka \
        --outdir $OUT/$base \
        --prefix $base \
        --genus Campylobacter \
        --usegenus \
        --cpus 60 \
        "$f" \
        > $LOG/prokka_$base.log 2>&1
done
