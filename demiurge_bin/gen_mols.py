import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ANSI color codes
COLORS = [
    '\033[38;5;46m',   # Green
    '\033[38;5;196m',  # Red
    '\033[38;5;214m'   # Orange
]
RESET = '\033[0m'

def standardize_smiles(smiles: str) -> str:
    """
    Converts the input SMILES string to an RDKit canonical SMILES.
    Raises ValueError if the SMILES cannot be parsed by RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)

def print_progress(current: int, total: int) -> None:
    """
    Prints a dynamic progress bar with color to indicate progress.

    :param current: The current file being processed.
    :param total:   The total number of files to process.
    """
    bar_length = 25  # Length of the progress bar
    filled_length = int(bar_length * (current / total))

    # Build the progress bar with colored blocks
    color_cycle = COLORS[0]  # Green color for progress
    bar = color_cycle + '█' * filled_length + '-' * (bar_length - filled_length) + RESET

    percent = int(100 * current / total)
    sys.stdout.write(f"\rProgress: |{bar}| {current}/{total} ({percent}%)")
    sys.stdout.flush()

    if current == total:
        print('')  # Add a newline after the last update

def generate_mol_files(csv_path: str) -> str:
    """
    Generates .mol files from SMILES strings provided in a CSV file.
    Fails are logged in 'mol_creation_error.log' and no .mol file is created for them.

    :param csv_path: Path to the input CSV file containing 'MOLECULE_NAME' and 'SMILES' columns.
    :return: The path to the directory where .mol files are stored.
    """
    mols_directory = os.path.join(os.getcwd(), "mols")
    error_log_path = os.path.join(os.getcwd(), "mol_creation_error.log")

    # Create the 'mols' directory if it does not exist
    if not os.path.exists(mols_directory):
        os.makedirs(mols_directory)
        print(f"\nCreated directory: {COLORS[2]}{mols_directory}{RESET}")

    try:
        data = pd.read_csv(csv_path)
        if 'MOLECULE_NAME' not in data.columns or 'SMILES' not in data.columns:
            raise ValueError(f"{COLORS[1]}CSV must contain 'MOLECULE_NAME' and 'SMILES' columns.{RESET}")

        # Remove duplicates in the 'MOLECULE_NAME' column
        initial_count = len(data)
        data = data.drop_duplicates(subset='MOLECULE_NAME', keep='first')
        duplicates_count = initial_count - len(data)
        if duplicates_count > 0:
            print(f"{COLORS[1]}\nRemoved {duplicates_count} duplicate entries in 'MOLECULE_NAME'.{RESET}")

        # Counters
        file_count = 0        # Number of successful .mol files
        error_file_count = 0  # Number of failed .mol files
        errors = []           # List to store error details

        print("\nGenerating *.mol files ...\n")

        total_files = len(data)
        last_update = 0

        for index, row_data in enumerate(data.iterrows(), 1):
            name = row_data[1]['MOLECULE_NAME']
            raw_smiles = row_data[1]['SMILES']

            try:
                # 1) Standardize the SMILES
                smiles = standardize_smiles(raw_smiles)

                # 2) Create the molecule from the standardized SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Could not parse standardized SMILES: {smiles}")

                # Add hydrogen atoms
                mol = Chem.AddHs(mol)

                # 3) Generate 3D coordinates (for stereochemistry and other computations)
                AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
                # Optionally, you could also do a force-field optimization
                # AllChem.MMFFOptimizeMolecule(mol)

                # 4) Generate 2D coordinates (overwrites conformer with 2D layout, 
                #    but keeps the stereochemical info in the graph)
                AllChem.Compute2DCoords(mol)

                # Remove partial double-bond stereo for any N=C bonds
                for bond in mol.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        a1 = bond.GetBeginAtom()
                        a2 = bond.GetEndAtom()
                        # If any of these atoms is nitrogen, clear stereo info
                        if a1.GetAtomicNum() == 7 or a2.GetAtomicNum() == 7:
                            bond.SetStereo(Chem.BondStereo.STEREONONE)
                            # Removed: bond.SetStereoAtoms(-1, -1)

                # Write .mol file (V3000 format)
                mol_file_path = os.path.join(mols_directory, f"{name}.mol")
                with open(mol_file_path, 'w') as f:
                    f.write(Chem.MolToMolBlock(mol, forceV3000=True))

                file_count += 1

            except Exception as e:
                error_file_count += 1
                errors.append(
                    f"Molecule: {name}, SMILES: {raw_smiles}\nError: {str(e)}\n"
                )

            # Update progress bar
            progress = (index / total_files) * 100
            if index != total_files and progress - last_update >= 1:
                print_progress(index, total_files)
                last_update = progress

        # Final update to 100%
        print_progress(total_files, total_files)

        # If any errors occurred, write them to the log file
        if error_file_count > 0:
            with open(error_log_path, 'w', encoding='utf-8') as logf:
                logf.write("==== MOL CREATION ERRORS ====\n\n")
                for err in errors:
                    logf.write(err + "\n")

        print(f"\n{COLORS[0]}Generated {file_count} .mol files in '{mols_directory}'.{RESET}")
        print(f"{COLORS[0]}Failed to generate {error_file_count} .mol files.{RESET}")
        if error_file_count > 0:
            print(f"{COLORS[1]}Details are in the log file: '{error_log_path}'{RESET}")

    except FileNotFoundError:
        print(f"{COLORS[1]}File not found: {csv_path}{RESET}")
    except pd.errors.EmptyDataError:
        print(f"{COLORS[1]}File is empty: {csv_path}{RESET}")
    except pd.errors.ParserError:
        print(f"{COLORS[1]}Error parsing CSV file: {csv_path}{RESET}")
    except Exception as e:
        print(f"{COLORS[1]}An unexpected error occurred: {e}{RESET}")

    return mols_directory
