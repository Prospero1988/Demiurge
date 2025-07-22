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
import re
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
from demiurge_bin.fp_generator import fp_generator

# Import sys and define the Tee class
import sys

def strip_ansi_codes(s):
    ansi_escape = re.compile(r'''
        \x1B # ESC
        (?:   # 7-bit C1 Fe (various sequences)
            [@-Z\\-_]
        | # or CSI [ - ].
            \[
            [0-?]* # Optional parameters.
            [-/]* # Optional intermediate bytes.
            [@-~] # Final byte
        )
    ''', re.VERBOSE)
    return ansi_escape.sub('', s)

# Import sys and define the Tee class
class Tee(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            # Check if the file is a terminal (stdout/stderr)
            if hasattr(f, 'isatty') and f.isatty():
                f.write(obj)
            else:
                # Remove ANSI codes before writing to file
                f.write(strip_ansi_codes(obj))
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()

def main():
    """
    Main function to execute the NMR data processing pipeline.
    """
    
    # ANSI color
    COLORS = ['\033[38;5;46m',    # Green
              '\033[38;5;196m',   # Red
              '\033[38;5;214m'    # Orange
             ]
    RESET = '\033[0m'
    
    # Open the log file in append mode
    log_file = open('demiurge.log', 'a', encoding='utf-8')

    # Redirect sys.stdout and sys.stderr to the Tee instance
    sys.stdout = Tee(sys.__stdout__, log_file)
    sys.stderr = Tee(sys.__stderr__, log_file)

    try:

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
            choices=["1H", "13C", "FP"],
            help="Type of NMR Predictor to use: '1H' or '13C'. You can also "
            "generate from RDKit FingerPrints not from NMR Spectra, "
            "use argument: 'FP'"
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
    
        # Clear the console and display the ASCII art logo
        subprocess.call('cls' if os.name == 'nt' else 'clear', shell=True)
        
        print('')
        ascii_art_demiurge = text2art("DEMIURGE")
        predictor = args.predictor
        ascii_art_predictor = text2art(predictor)
        art_width = len(ascii_art_demiurge.split('\n')[0])
        centered_predictor_lines = [line.center(art_width) for line in ascii_art_predictor.split('\n')]
        final_art = f"{ascii_art_demiurge}\n" + "\n".join(centered_predictor_lines)
        print(final_art)                   

        # Step 1: Verify the CSV input file and correct any issues
        verified_csv_path = verify_csv(args.csv_path)

        if args.predictor == 'FP':
            
            # Step 2 - 4: Generate FingerPrint files
            processed_dir = fp_generator(verified_csv_path)
            
        else:

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
                f"Script executed with the {COLORS[2]}--clean {RESET}option. All temporary files "
                f"and folders will be removed:\n"
            )
                
            temp_data_fp = [processed_dir, merged_dir]
            
            # Assign temp_data based on predictor
            if args.predictor != "FP":
                temp_data_fp.extend([csv_output_folder, mol_directory])
            temp_data = temp_data_fp

            if os.path.exists(verified_csv_path):
                os.remove(verified_csv_path)
                print(f"The file {COLORS[2]}'{verified_csv_path}'{RESET} has been deleted.")
            else:
                print(f"{COLORS[1]}The file '{verified_csv_path}' does not exist.{RESET}")
            
            for folder in temp_data:
                if os.path.exists(folder):
                    shutil.rmtree(folder)
                    print(f"Temporary folder {COLORS[2]}'{folder}'{RESET} has been deleted.")
                else:
                    print(f"{COLORS[1]}Folder '{folder}' does not exist.{RESET}")
    
    finally:
        # Restore original sys.stdout and sys.stderr
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        
        # Close the log file
        log_file.close()
    
    print(f"\n{COLORS[0]}END OF THE SCRIPT\n{RESET}")

if __name__ == "__main__":
    main()
