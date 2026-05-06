#!/bin/bash

# Function to display help message
display_help() {
    echo "Usage: $0 [options]"
    echo
    echo "This script runs MLST on all genome files in the specified folder and extracts"
    echo "only the results line for each genome, saving them to the specified output file."
    echo
    echo "Options:"
    echo "  -i, --input <input_folder>   Path to the folder containing genome files (.fna format)"
    echo "  -o, --output <output_file>   Name of the file to save the MLST results"
    echo "  -h, --help                   Display this help message and exit"
    echo
    echo "Example:"
    echo "  $0 -i Africa_genomes -o mlst_results.txt"
    echo "  $0 --input Africa_genomes --output mlst_results.txt"
    echo
    echo "Dependencies:"
    echo "  - 'mlst' tool must be installed and available in the PATH."
    echo "  - The folder must contain genome files in .fna format."
    echo
    exit 0
}

# Initialize variables for input folder and output file
input_folder=""
output_file=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input_folder="$2"
            shift 2
            ;;
        -o|--output)
            output_file="$2"
            shift 2
            ;;
        -h|--help)
            display_help
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage instructions."
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$input_folder" || -z "$output_file" ]]; then
    echo "Error: Missing required arguments."
    echo "Use -h or --help for usage instructions."
    exit 1
fi

# Ensure output file is empty before writing new results
> "$output_file"

# Loop through all genome files in the input folder
for genome in "$input_folder"/*.fna; do
    # Extract the base name of the genome file (without path and extension)
    genome_base=$(basename "$genome")

    # Run MLST on the current genome
    mlst_output=$(mlst --nopath "$genome")

    # Extract the line containing the MLST results for this genome (starting with the base name)
    echo "$mlst_output" | grep -E "^${genome_base}" >> "$output_file"
done

echo "MLST results saved to $output_file"
