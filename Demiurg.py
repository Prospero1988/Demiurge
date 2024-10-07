# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 20:09:08 2024

@author: arkadiusz.leniak@gmal.com
"""
# import of general libraries
import argparse
import shutil
import os

# import executable modules for this script
from demiurge_bin.gen_mols import generate_mol_files
from demiurge_bin.predictor import run_java_batch_processor
from demiurge_bin.bucket import bucket
from demiurge_bin.merger import merger
from demiurge_bin.labeler import labeler
from demiurge_bin.custom_header import custom_header

def main():
    # Define the argument parser
    parser = argparse.ArgumentParser(description="Complex software to create NMR-based ML inputs from SMILES")
    parser.add_argument("--csv_path", type=str, required=True, help="Path to the input CSV file with SMILES.")
    parser.add_argument("--predictor", type=str, required=True, help="Choose 1H or 13C for NMR Predictor.")
    parser.add_argument("--label_column", type=int, required=True, help="Specify column number (as integer) in input file in which labels are present.")
    parser.add_argument("--clean", action='store_true', help="If used, script deletes all intermediate temp files.")

    # Parse the arguments
    args = parser.parse_args()

    # Generate *.mol files from SMILES
    mol_directory = generate_mol_files(args.csv_path)

    # Predict NMR spectra and save as *.mol files
    csv_output_folder = run_java_batch_processor(mol_directory, args.predictor)

    # Bucketing process - generating pseudo NMR Spectra
    processed_dir = bucket(csv_output_folder, args.predictor)

    # Merging Spectra in CSV format to one matrix file.
    output_path, merged_dir = merger(processed_dir, args.csv_path)

    # Insert Labels to merged files.
    labeled = labeler(args.csv_path, output_path, args.label_column, merged_dir)

    # Creating Custom Headers
    custom_header(labeled, args.csv_path)

    # Cleaning and removing temporary dirs and data if --clean flag is set
    if args.clean:
        print("Script executed with option --clean. All temp files and folders will be removed:\n")
        temp_data = [mol_directory, csv_output_folder, processed_dir, merged_dir]
        for folder in temp_data:
            if os.path.exists(folder):
                shutil.rmtree(folder)
                print(f"Temp folder '{folder}' has been deleted.")
            else:
                print(f"Folder '{folder}' doesn't exist.")

if __name__ == "__main__":
    main()

print("\nEND OF THE SCRIPT\n")