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
    """Runs a shell command and optionally logs its output."""
    try:
        result = subprocess.run(
            command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if log_file:
            with open(log_file, "a") as log:
                log.write(result.stdout + "\n" + result.stderr + "\n")
        if result.returncode != 0:
            print_warning(f"Command failed: {command}\n{result.stderr}")
    except Exception as e:
        print_warning(f"Error executing command: {command}\n{e}")

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

def annotate_with_prokka(genome_dir, output_dir, threads, kingdom):
    """Run Prokka annotation on genomes."""
    print_info("Starting genome annotation with Prokka...")
    prokka_output_dir = os.path.join(output_dir, "prokka_results")
    os.makedirs(prokka_output_dir, exist_ok=True)
    
    for genome in os.listdir(genome_dir):
        if genome.endswith(".fna"):
            base_name = genome.split(".")[0]
            genome_path = os.path.join(genome_dir, genome)
            genome_output = os.path.join(prokka_output_dir, base_name)
            
            print_info(f"Annotating {base_name} with Prokka...")
            prokka_cmd = f"prokka --cpus {threads} --kingdom {kingdom} --outdir {genome_output} --force --norrna --notrna --addgenes  {genome_path}"
            log_file = os.path.join(output_dir, f"{base_name}_prokka.log")
            run_command(prokka_cmd, log_file)

    print_info("Genome annotation with Prokka completed.")

def create_sample_sheet(prokka_results_dir, output_dir):
    """Create a sample sheet for Abricate."""
    print_info("Creating sample sheet for AbritAMR...")
    sample_sheet = os.path.join(output_dir, "sample_sheet.txt")
    
    with open(sample_sheet, "w") as sheet:
        for genome_folder in os.listdir(prokka_results_dir):
            assembly_file = os.path.join(prokka_results_dir, genome_folder, "PROKKA.fna")
            if os.path.isfile(assembly_file):
                sheet.write(f"{genome_folder}\t{assembly_file}\n")

    print_info(f"Sample sheet created at {sample_sheet}")
    return sample_sheet

def run_abritamr(sample_sheet, output_dir, threads, species):
    """Run AbritAMR for AMR prediction."""
    print_info("Starting AMR prediction with AbritAMR...")
    abritamr_output_dir = os.path.join(output_dir, "abritamr_results")
    os.makedirs(abritamr_output_dir, exist_ok=True)
    
    abritamr_cmd = f"abritamr run -j {threads} --species {species} -c {sample_sheet}"
    log_file = os.path.join(output_dir, "abritamr.log")
    run_command(abritamr_cmd, log_file)

    print_info("Moving AbritAMR results to output folder...")
    for result_file in ["summary_matches.txt", "summary_partials.txt", "summary_virulence.txt"]:
        if os.path.isfile(result_file):
            os.rename(result_file, os.path.join(abritamr_output_dir, result_file))

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
        "--species", "-s",
        required=True,
        choices=[
            "Campylobacter", "Escherichia", "Klebsiella_pneumoniae", "Salmonella", 
            "Staphylococcus_aureus", "Vibrio_cholerae"
        ],
        help="Species for AbritAMR AMR prediction. Required for point mutation analysis."
    )
    args = parser.parse_args()

    # Pipeline Execution
    start_time = datetime.now()
    print_info("Pipeline started.")
    print_info(f"Input directory: {args.input_dir}")
    print_info(f"Output directory: {args.output_dir}")
    print_info(f"Using {args.threads} threads.")
    print_info(f"Kingdom for Prokka: {args.kingdom}")
    print_info(f"Species for AbritAMR: {args.species}")

    # Step 1: Decompress .gz Files
    decompress_gz_files(args.input_dir)

    # Step 2: Run Prokka
    annotate_with_prokka(args.input_dir, args.output_dir, args.threads, args.kingdom)

    # Step 3: Create Sample Sheet
    prokka_results_dir = os.path.join(args.output_dir, "prokka_results")
    sample_sheet = create_sample_sheet(prokka_results_dir, args.output_dir)

    # Step 4: Run AbritAMR
    run_abritamr(sample_sheet, args.output_dir, args.threads, args.species)

    # Completion
    end_time = datetime.now()
    print_info(f"Pipeline completed in {end_time - start_time}.")

if __name__ == "__main__":
    main()
