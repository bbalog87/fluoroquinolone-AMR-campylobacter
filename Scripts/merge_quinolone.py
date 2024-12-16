import os
import argparse
import pandas as pd

def merge_quinolone(input_dir, output_file):
    """
    Merge all quinolone result files in the input directory into a single master file.
    Ensures all unique columns are included in the master file, with missing values filled as NaN.
    
    Args:
        input_dir (str): Path to the directory containing input files.
        output_file (str): Path to the output file.
    """
    # Validate input directory
    if not os.path.isdir(input_dir):
        print(f"Error: Directory '{input_dir}' not found.")
        return
    
    # List all .txt files in the input directory
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.txt')]
    if not files:
        print(f"No .txt files found in directory '{input_dir}'.")
        return
    
    # Initialize an empty DataFrame for merging
    master_df = pd.DataFrame()

    for file in files:
        print(f"Processing file: {file}")
        try:
            # Read the current file
            current_df = pd.read_csv(file, sep="\t")
            
            # Merge with the master DataFrame, aligning on columns
            master_df = pd.concat([master_df, current_df], ignore_index=True, sort=False)
        except Exception as e:
            print(f"Error processing file '{file}': {e}")
    
    # Save the merged DataFrame to the output file
    print(f"Saving merged data to: {output_file}")
    master_df.to_csv(output_file, sep="\t", index=False)
    print("Merge complete.")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Merge all quinolone result files in a directory into a single master file."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the directory containing input quinolone result files."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path to the output file for the merged results."
    )
    
    args = parser.parse_args()

    # Call the merge function
    merge_quinolone(args.input, args.output)
