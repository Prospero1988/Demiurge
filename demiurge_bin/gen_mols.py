#gen_mols.py

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_mol_files(csv_path):
    """
    Function to generate .mol files from SMILES strings provided in a CSV file.
    Attempts to fix invalid SMILES structures using RDKit functions.

    Parameters:
    - csv_path: Path to the input CSV file containing 'MOLECULE_NAME' and 'SMILES' columns.

    Returns:
    - mols_directory: Path to the directory where generated .mol files are stored.
    """
    
    # Define the path to the "mols" subdirectory relative to the main script directory
    mols_directory = os.path.join(os.getcwd(), "mols")
    
    # Check if the "mols" directory exists; if not, create it
    if not os.path.exists(mols_directory):
        os.makedirs(mols_directory)
        print(f"\nCreated directory: {mols_directory}")
    
    try:
        # Load the CSV file
        data = pd.read_csv(csv_path)
        
        # Check if the necessary columns exist in the CSV file
        if 'MOLECULE_NAME' not in data.columns or 'SMILES' not in data.columns:
            raise ValueError("CSV file must contain 'MOLECULE_NAME' and 'SMILES' columns.")
        
        file_count = 0
        
        print("\nGenerating *.mol files ...")
        
        # Iterate through each row of the CSV file
        for index, row in data.iterrows():
            try:
                name = row['MOLECULE_NAME']
                smiles = row['SMILES']
                
                # Generate the molecule from the SMILES string
                mol = Chem.MolFromSmiles(smiles)
                
                # If the molecule is invalid, try to sanitize and kekulize the molecule
                if mol is None:
                    print(f"Invalid SMILES string at row {index}: {smiles}. Trying to fix...")
                    mol = Chem.MolFromSmiles(smiles, sanitize=False)
                    if mol:
                        try:
                            # Try sanitizing the molecule
                            Chem.SanitizeMol(mol)
                            print(f"Successfully sanitized molecule at row {index}.")
                        except:
                            print(f"Sanitization failed for molecule at row {index}. Attempting kekulization...")
                            try:
                                Chem.Kekulize(mol, clearAromaticFlags=True)
                                print(f"Successfully kekulized molecule at row {index}.")
                            except Exception as kek_error:
                                print(f"Kekulization failed for molecule at row {index}: {kek_error}")
                                continue

                if mol is None:
                    raise ValueError(f"Failed to generate molecule from SMILES: {smiles}")

                # Add hydrogen atoms to the molecule
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates for the molecule
                AllChem.EmbedMolecule(mol, randomSeed=42)  # EmbedMolecule can fail randomly, so use a seed
                
                # Optimize the 3D geometry
                AllChem.UFFOptimizeMolecule(mol)
                
                # Save the molecule to a .mol file in the "mols" subdirectory
                mol_file_path = os.path.join(mols_directory, f"{name}.mol")
                with open(mol_file_path, 'w') as f:
                    f.write(Chem.MolToMolBlock(mol))
                
                file_count += 1
                
            except Exception as e:
                print(f"Error processing row {index}: {e}")
        
        print(f"\nGenerated {file_count} .mol files in '{mols_directory}'")
    
    except FileNotFoundError:
        print(f"File not found: {csv_path}")
    except pd.errors.EmptyDataError:
        print(f"File is empty: {csv_path}")
    except pd.errors.ParserError:
        print(f"Error parsing CSV file: {csv_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
    return mols_directory

# Example usage:
# csv_file_path = "path_to_your_smiles_file.csv"
# generate_mol_files(csv_file_path)
