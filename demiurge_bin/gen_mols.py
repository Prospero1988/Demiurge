#gen_mols.py

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_mol_files(csv_path):
    """
    Function to generate .mol files from SMILES strings provided in a CSV file.

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
        
        # Iterate through each row of the CSV file
        for index, row in data.iterrows():
            try:
                name = row['MOLECULE_NAME']
                smiles = row['SMILES']
                
                # Generate the molecule from the SMILES string
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Invalid SMILES string: {smiles}")
                
                # Add hydrogen atoms to the molecule
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates for the molecule
                AllChem.EmbedMolecule(mol)
                
                # Optionally, optimize the molecule in 2D
                AllChem.Compute2DCoords(mol)
                
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
