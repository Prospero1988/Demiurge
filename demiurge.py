#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script Name: demiurge.py
Description: A comprehensive software pipeline to create NMR-based
machine learning inputs from SMILES strings.
"""

import argparse
import shutil
import os
import subprocess
from art import text2art

# Import custom modules required for the script
from demiurge_bin.csv_checker import verify_csv
from demiurge_bin.gen_mols import generate_mol_files
from demiurge_bin.predictor import run_java_batch_processor
from demiurge_bin.bucket import bucket
from demiurge_bin.merger import merger
from demiurge_bin.labeler import labeler
from demiurge_bin.custom_header import custom_header


def main():
    """
    Main function to execute the NMR data processing pipeline.
    """

    # Clear the console and display the ASCII art logo
    subprocess.call('cls' if os.name == 'nt' else 'clear', shell=True)
    print('')
    ascii_art = text2art("DEMIURGE")
    print(ascii_art)

    # Define the argument parser for command-line options
    parser = argparse.ArgumentParser(
        description="A comprehensive software pipeline to create "
                    "NMR-based machine learning inputs from SMILES strings."
    )
    parser.add_argument(
        "--csv_path",
        type=str,
        required=True,
        help="Path to the input CSV file containing SMILES strings."
    )
    parser.add_argument(
        "--predictor",
        type=str,
        required=True,
        choices=["1H", "13C"],
        help="Type of NMR Predictor to use: '1H' or '13C'."
    )
    parser.add_argument(
        "--label_column",
        type=int,
        required=True,
        help="Column number (as an integer) in the input CSV file "
             "where labels are present."
    )
    parser.add_argument(
        "--clean",
        action='store_true',
        help="If set, the script deletes all intermediate temporary "
             "files after execution."
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Step 1: Verify the CSV input file and correct any issues
    verified_csv_path = verify_csv(args.csv_path)

    # Step 2: Generate .mol files from SMILES strings
    mol_directory = generate_mol_files(verified_csv_path)

    # Step 3: Predict NMR spectra and save results as .csv files
    csv_output_folder = run_java_batch_processor(mol_directory, args.predictor)

    # Step 4: Perform bucketing to generate pseudo NMR spectra
    processed_dir = bucket(csv_output_folder, args.predictor)

    # Step 5: Merge spectra in CSV format into one matrix file
    output_path, merged_dir = merger(processed_dir, verified_csv_path)

    # Step 6: Insert labels into the merged files
    labeled = labeler(
        verified_csv_path, output_path, args.label_column, merged_dir
    )

    # Step 7: Create custom headers for the final dataset
    custom_header(labeled, verified_csv_path, args.predictor)

    # Optional: Clean up temporary dirs and data if the --clean flag is set
    if args.clean:
        print(
            "Script executed with the --clean option. All temporary files "
            "and folders will be removed:\n"
        )
        temp_data = [
            mol_directory, csv_output_folder, processed_dir, merged_dir
        ]
        for folder in temp_data:
            if os.path.exists(folder):
                shutil.rmtree(folder)
                print(f"Temporary folder '{folder}' has been deleted.")
            else:
                print(f"Folder '{folder}' does not exist.")


if __name__ == "__main__":
    main()

print("\nEND OF THE SCRIPT\n")
