#!/usr/bin/env bash
#
# Copy .fna genome files whose filenames begin with accessions listed
# in accessions.txt. Files are copied from genomes/ to selected_genomes/.
#
# Assumptions:
# - Each accession appears at the start of the filename.
# - Filenames look like: GCA_XXXXXXX.Y_genomic.fna
# - Only .fna files should be copied.

set -euo pipefail

ACC_FILE="extended_failed.txt"
SOURCE_DIR="all_genomes"
TARGET_DIR="extended_failed_genomes"

mkdir -p "$TARGET_DIR"

# Read accessions line by line
while IFS= read -r acc; do
    # Skip empty lines
    [[ -z "$acc" ]] && continue

    # Pattern: accession at start, anything after, ending in .fna
    for file in "$SOURCE_DIR"/"$acc"*.fna; do
        # If no file matches, skip
        [[ -e "$file" ]] || continue

        cp -p "$file" "$TARGET_DIR"/
        echo "Copied: $(basename "$file")"
    done
done < "$ACC_FILE"
