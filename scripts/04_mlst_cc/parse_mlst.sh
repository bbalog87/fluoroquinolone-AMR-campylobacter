#!/bin/bash

# Usage: ./parse_mlst.sh <input_file> <output_file>
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

# Input and output files
input_file="$1"
output_file="$2"

# Temporary file for storing unique gene names
temp_genes_file=$(mktemp)

# Step 1: Collect all unique genes
grep -oP '\b\w+_\w+\(\d+\)|\b\w+\(\d+\)' "$input_file" | awk -F'(' '{print $1}' | sort -u > "$temp_genes_file"

# Convert the list of genes into a comma-separated header
genes=$(tr '\n' ',' < "$temp_genes_file" | sed 's/,$//')

# Step 2: Write the header to the output file
echo "Isolate,Species,ST,$genes" > "$output_file"

# Step 3: Process each line of the input file
while IFS= read -r line; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Extract isolate name, species, and ST
    isolate=$(echo "$line" | awk '{print $1}' | sed 's/\.fna$//')
    species=$(echo "$line" | awk '{print $2}')
    st=$(echo "$line" | awk '{print $3}')

    # Initialize an associative array for gene alleles
    declare -A gene_alleles

    # Extract all gene-allele pairs
    while read -r gene; do
        allele=$(echo "$line" | grep -oP "$gene\(\K[^\)]*" || echo "")
        gene_alleles["$gene"]="$allele"
    done < "$temp_genes_file"

    # Construct the row for this isolate
    row="$isolate,$species,$st"
    while read -r gene; do
        row="$row,${gene_alleles[$gene]:-}"  # Use empty value if gene is absent
    done < "$temp_genes_file"

    # Write the row to the output file
    echo "$row" >> "$output_file"

    # Clear the associative array
    unset gene_alleles
done < "$input_file"

# Cleanup
rm -f "$temp_genes_file"

echo "Parsing complete. Results saved to $output_file"
