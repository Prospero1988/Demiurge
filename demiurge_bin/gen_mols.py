# gen_mnols.py

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_mol_files(csv_path):
    # Ustal ścieżkę do podkatalogu "mols" względem głównego skryptu
    mols_directory = os.path.join(os.getcwd(), "mols")
    
    # Sprawdź, czy katalog "mols" istnieje, jeśli nie, to go utwórz
    if not os.path.exists(mols_directory):
        os.makedirs(mols_directory)
        print(f"\nCreated directory: {mols_directory}")
    
    try:
        # Wczytaj plik CSV
        data = pd.read_csv(csv_path)
        
        # Sprawdź, czy wymagane kolumny istnieją
        if 'MOLECULE_NAME' not in data.columns or 'SMILES' not in data.columns:
            raise ValueError("CSV file must contain 'MOLECULE_NAME' and 'SMILES' columns.")
        
        file_count = 0
        
        for index, row in data.iterrows():
            try:
                name = row['MOLECULE_NAME']
                smiles = row['SMILES']
                
                # Wygeneruj cząsteczkę z SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Invalid SMILES string: {smiles}")
                
                # Dodaj atomy wodoru
                mol = Chem.AddHs(mol)
                
                # Wygeneruj współrzędne 3D
                AllChem.EmbedMolecule(mol)
                
                # Optymalizuj cząsteczkę w 2D
                #AllChem.Compute2DCoords(mol)
                
                # Zapisz cząsteczkę do pliku .mol w podkatalogu "mols"
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