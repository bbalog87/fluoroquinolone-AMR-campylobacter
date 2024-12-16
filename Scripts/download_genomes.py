import os
import sys
from Bio import Entrez
import argparse

# Set your email (required by NCBI to identify users)
Entrez.email = "your_email@example.com"

def download_genome(accession, output_dir):
    """Download genome sequence from NCBI by genome assembly accession number.

    Args:
        accession (str): The genome assembly accession number.
        output_dir (str): The directory to save downloaded genome files.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Search for the genome in the Assembly database
        search_handle = Entrez.esearch(db="assembly", term=accession, retmode="xml")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        # Get the assembly UID
        if search_results['IdList']:
            assembly_uid = search_results['IdList'][0]
            
            # Fetch the assembly summary
            summary_handle = Entrez.esummary(db="assembly", id=assembly_uid, retmode="xml")
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            
            # Retrieve the FTP link for the genome sequence
            ftp_link = summary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
            if ftp_link:
                fasta_link = f"{ftp_link}/{ftp_link.split('/')[-1]}_genomic.fna.gz"
                print(f"Downloading genome for {accession} from {fasta_link}")
                
                # Download the genome sequence
                genome_file = os.path.join(output_dir, f"{accession}.fna.gz")
                os.system(f"wget -O {genome_file} {fasta_link}")
                print(f"Downloaded: {genome_file}")
            else:
                print(f"No FTP link found for {accession}")
        else:
            print(f"No results found for accession {accession}")
    except Exception as e:
        print(f"Error downloading genome for {accession}: {e}")

def main(input_file, output_dir):
    """Main function to download genomes from a list of accession numbers."""
    # Read the list of accession numbers from the input file
    with open(input_file, 'r') as file:
        accessions = [line.strip() for line in file if line.strip()]
    
    # Download each genome based on accession number
    for accession in accessions:
        download_genome(accession, output_dir)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Download genome sequences from NCBI based on assembly accession numbers.",
        epilog="Example: python download_genomes.py -i accessions.txt -o genomes"
    )
    
    # Add arguments for input file and output directory
    parser.add_argument("-i", "--input_file", required=True, help="File containing genome assembly accession numbers (one per line).")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save the downloaded genome files.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run the main function
    main(args.input_file, args.output_dir)
