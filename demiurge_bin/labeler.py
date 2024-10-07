# labeler.py

import pandas as pd
import os

def labeler(csv_path, output_path, label_column, merged_dir):
    # Wczytanie pliku CSV z etykietami i SMILES
    labels_csv = pd.read_csv(csv_path)

    # Wczytanie scalonego pliku CSV z danymi NMR
    merged_nmr = pd.read_csv(output_path)
    
    # Dodanie nowej kolumny 'LABEL' z początkową wartością 0
    merged_nmr['LABEL'] = 0

    # Przeprowadzenie dopasowania na podstawie nazw plików
    labeled = pd.merge(merged_nmr, labels_csv[['MOLECULE_NAME', label_column]], 
                      left_on='filename', right_on='MOLECULE_NAME', how='left')

    # Przypisanie wartości z kolumny `label_column` do `LABEL`
    labeled['LABEL'] = labeled[label_column]

    # Usunięcie zbędnej kolumny `MOLECULE_NAME`, jeśli chcesz
    labeled.drop(columns=['MOLECULE_NAME'], inplace=True)

    file_name = os.path.basename(csv_path).split('.')[0] + '_labeled.csv'
    labeled.to_csv(file_name, index=False)
    
    return labeled

