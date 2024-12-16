import argparse
import pandas as pd

def extract_quinolone(input_file, output_file):
    """
    Extract quinolone-related resistance information from AbritAMR results.
    
    Parameters:
    - input_file: Path to the input AbritAMR results file (tab-delimited).
    - output_file: Path to the output file for extracted data.
    """
    try:
        # Load the input file as a tab-delimited DataFrame
        print(f"Loading data from: {input_file}")
        data = pd.read_csv(input_file, sep="\t")

        # Identify columns containing "Quinolone" (case-insensitive)
        quinolone_columns = [col for col in data.columns if "quinolone" in col.lower()]

        # Validate that we found relevant columns
        if not quinolone_columns:
            print("Error: No columns containing 'Quinolone' found in the input file.")
            return

        # Select relevant columns (Isolate + Quinolone-related columns)
        selected_columns = ["Isolate"] + quinolone_columns
        extracted_data = data[selected_columns]

        # Save the extracted data to the output file
        print(f"Saving extracted data to: {output_file}")
        extracted_data.to_csv(output_file, sep="\t", index=False)
        print("Extraction complete.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description="""
        Extract Quinolone-related resistance information from AbritAMR result files.
        
        This script identifies columns containing the word 'Quinolone' (case-insensitive) and 
        outputs a filtered file with the selected rows and relevant columns for downstream analysis.
        
        The script assumes the first row of the input file is the header row.
        """,
        epilog="""
        Example Usage:
          python extract_quinolone.py -i abritamr_results.txt -o quinolone_results.txt

        Notes:
          - Empty cells in the 'Quinolone' columns indicate no resistance.
          - Non-empty cells (containing gene names) indicate quinolone resistance.
          - The script saves the output as a tab-delimited file.
        """
    )
    parser.add_argument(
        "-i", "--input", required=True, 
        help="Path to the input AbritAMR result file (tab-delimited)."
    )
    parser.add_argument(
        "-o", "--output", required=True, 
        help="Path to the output file for extracted data."
    )

    # Parse arguments
    args = parser.parse_args()

    # Extract quinolone resistance information
    extract_quinolone(args.input, args.output)
