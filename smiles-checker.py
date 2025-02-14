import argparse
import pandas as pd
from rdkit import Chem

def check_smiles_validity(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    return True, "Valid SMILES"

def process_csv(input_csv):
    df = pd.read_csv(input_csv)
    if len(df.columns) < 2:
        print("Error: The input CSV must have at least two columns.")
        return
    
    errors = []
    for idx, row in df.iterrows():
        smiles = str(row.iloc[1]).strip()
        is_valid, message = check_smiles_validity(smiles)
        if not is_valid or "stereochemistry" in message:
            errors.append((idx + 1, smiles, message))
    
    if errors:
        print("Errors found in the following rows:")
        for row, smiles, msg in errors:
            print(f"Row {row}: {smiles} -> {msg}")
    else:
        print("All SMILES structures are valid.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check validity of SMILES structures in a CSV file.")
    parser.add_argument("input_csv", help="Path to the input CSV file.")
    args = parser.parse_args()
    
    process_csv(args.input_csv)
