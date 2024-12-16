#!/usr/bin/env python3

import os
import argparse
import subprocess
import gzip
import shutil
from datetime import datetime


# Helper Functions
def print_info(message):
    print(f"[INFO] {message}")


def print_warning(message):
    print(f"[WARNING] {message}")


def run_command(command, log_file=None):
    """Run a shell command and optionally log its output."""
    try:
        result = subprocess.run(
            command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if log_file:
            with open(log_file, "a") as log:
                log.write(result.stdout + "\n" + result.stderr + "\n")
        if result.returncode != 0:
            print_warning(f"Command failed: {command}\n{result.stderr}")
        return result.stdout
    except Exception as e:
        print_warning(f"Error executing command: {command}\n{e}")
        raise


def decompress_gz_files(genome_dir):
    """Decompress .gz files in the specified directory."""
    print_info("Checking for compressed .gz files...")
    for file in os.listdir(genome_dir):
        if file.endswith(".gz"):
            gz_file_path = os.path.join(genome_dir, file)
            decompressed_file_path = os.path.join(genome_dir, os.path.splitext(file)[0])
            
            print_info(f"Decompressing {file}...")
            with gzip.open(gz_file_path, 'rb') as gz_file:
                with open(decompressed_file_path, 'wb') as decompressed_file:
                    shutil.copyfileobj(gz_file, decompressed_file)
            os.remove(gz_file_path)  # Remove the .gz file after decompression
            print_info(f"Decompressed {file} to {decompressed_file_path}")
    print_info("Decompression of .gz files completed.")


def annotate_with_prokka(input_dir, output_dir, threads, kingdom):
    """Run Prokka annotation on genomes."""
    print_info("Starting genome annotation with Prokka...")
    prokka_output_dir = os.path.join(output_dir, "prokka_results")
    os.makedirs(prokka_output_dir, exist_ok=True)

    for genome_file in os.listdir(input_dir):
        genome_path = os.path.join(input_dir, genome_file)
        if os.path.isfile(genome_path) and genome_file.endswith(".fna"):
            base_name = os.path.splitext(genome_file)[0]
            genome_output = os.path.join(prokka_output_dir, base_name)
            
            print_info(f"Annotating {base_name} with Prokka...")
            prokka_cmd = (
                f"prokka --cpus {threads} --kingdom {kingdom} --outdir {genome_output} "
                f"--force --norrna --notrna {genome_path}"
            )
            log_file = os.path.join(output_dir, f"{base_name}_prokka.log")
            run_command(prokka_cmd, log_file)

            # Rename Prokka output files with the base name for clarity
            for prokka_file in os.listdir(genome_output):
                old_path = os.path.join(genome_output, prokka_file)
                new_path = os.path.join(genome_output, f"{base_name}_{prokka_file}")
                os.rename(old_path, new_path)

    print_info("Genome annotation with Prokka completed.")


def create_sample_sheet(input_dir, output_dir):
    """Create a sample sheet for AbritAMR."""
    print_info("Creating sample sheet for AbritAMR...")
    sample_sheet = os.path.join(output_dir, "sample_sheet.txt")
    entries = 0

    with open(sample_sheet, "w") as sheet:
        for genome_file in os.listdir(input_dir):
            genome_path = os.path.join(input_dir, genome_file)
            if os.path.isfile(genome_path) and genome_file.endswith(".fna"):
                abs_path = os.path.abspath(genome_path)
                base_name = os.path.splitext(genome_file)[0]
                sheet.write(f"{base_name}\t{abs_path}\n")
                entries += 1
            else:
                print_warning(f"Skipping {genome_file}: Not a valid .fna file.")

    if entries == 0:
        print_warning("No valid entries found. Sample sheet is empty.")
    else:
        print_info(f"Sample sheet created with {entries} entries at {sample_sheet}.")
    return sample_sheet


def run_abritamr(sample_sheet, output_dir, threads, species):
    """Run AbritAMR for AMR prediction and clean up base-name directories."""
    print_info("Starting AMR prediction with AbritAMR...")
    abritamr_output_dir = os.path.join(output_dir, "abritamr_results")
    os.makedirs(abritamr_output_dir, exist_ok=True)

    # Run AbritAMR
    abritamr_cmd = f"abritamr run -j {threads} --species {species} -c {sample_sheet}"
    log_file = os.path.join(output_dir, "abritamr.log")
    run_command(abritamr_cmd, log_file)

    # Move results to output directory
    for result_file in ["summary_matches.txt", "abritamr.txt", "summary_partials.txt", "summary_virulence.txt"]:
        if os.path.isfile(result_file):
            result_path = os.path.join(abritamr_output_dir, result_file)
            shutil.move(result_file, result_path)

    print_info("Attempting to clean up AbritAMR base-name directories...")

    # Identify and remove all directories corresponding to sample base names at the current level
    with open(sample_sheet, 'r') as sheet:
        sample_base_names = [line.split('\t')[0] for line in sheet if line.strip()]

    current_level_directories = [
        folder for folder in os.listdir(".")
        if os.path.isdir(folder) and folder in sample_base_names
    ]

    for folder in current_level_directories:
        try:
            shutil.rmtree(folder)
            print_info(f"Removed AbritAMR directory: {folder}")
        except Exception as e:
            print_warning(f"Failed to remove {folder}: {e}")

    print_info("Cleanup of base-name directories completed.")
    print_info("AMR prediction with AbritAMR completed.")



def main():
    # Argument Parsing
    parser = argparse.ArgumentParser(
        description="Pipeline for genome annotation and AMR prediction using Prokka and AbritAMR."
    )
    parser.add_argument(
        "--input_dir", "-i", required=True, help="Directory containing genome files in .fna or .fna.gz format."
    )
    parser.add_argument(
        "--output_dir", "-o", required=True, help="Directory where output files will be saved."
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=4, help="Number of threads to use. Default is 4."
    )
    parser.add_argument(
        "--kingdom", "-k", default="Bacteria", choices=["Bacteria", "Archaea", "Viruses", "Mitochondria", "Plasmids"],
        help="Kingdom for Prokka annotation. Default is 'Bacteria'."
    )
    parser.add_argument(
        "--species", "-s", required=True, help="Species for AbritAMR prediction (e.g., Campylobacter)."
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Pipeline Execution
    start_time = datetime.now()
    print_info("Pipeline started.")
    print_info(f"Input directory: {args.input_dir}")
    print_info(f"Output directory: {args.output_dir}")
    print_info(f"Using {args.threads} threads.")
    print_info(f"Kingdom for Prokka: {args.kingdom}")
    print_info(f"Species for AbritAMR: {args.species}")

    # Step 1: Decompress .gz files if necessary
    decompress_gz_files(args.input_dir)

    # Step 2: Run Prokka for annotation
#    annotate_with_prokka(args.input_dir, args.output_dir, args.threads, args.kingdom)

    # Step 3: Create Sample Sheet
    sample_sheet_path = create_sample_sheet(args.input_dir, args.output_dir)

    # Step 4: Run AbritAMR
    run_abritamr(sample_sheet_path, args.output_dir, args.threads, args.species)

    # Completion
    end_time = datetime.now()
    print_info(f"Pipeline completed in {end_time - start_time}.")

    print_info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
