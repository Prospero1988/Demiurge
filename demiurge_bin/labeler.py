# labeler.py

import pandas as pd
import os

def labeler(csv_path, output_path, label_column, merged_dir):
    try:
        # Wczytanie pliku CSV z etykietami i SMILES
        labels_csv = pd.read_csv(csv_path)
    
        # Wczytanie scalonego pliku CSV z danymi NMR
        merged_nmr = pd.read_csv(output_path)
        
        # Standaryzacja na pozycję 0 (zmiana label_column z 1-indeksowego na 0-indeksowe)
        label_column = label_column - 1
    
        # Pobranie nazwy kolumny na podstawie indeksu `label_column`
        label_column_name = labels_csv.columns[label_column]
    
        # Przeprowadzenie dopasowania na podstawie nazw plików (łączenie)
        labeled = pd.merge(merged_nmr, labels_csv[['MOLECULE_NAME', label_column_name]], 
                           left_on='filename', right_on='MOLECULE_NAME', how='left')
    
        # Przypisanie wartości z kolumny `label_column_name` do `LABEL`
        labeled['LABEL'] = labeled[label_column_name]
    
        # Usunięcie zbędnej kolumny `MOLECULE_NAME`
        labeled.drop(columns=['MOLECULE_NAME', label_column_name], inplace=True)
    
        # Przeniesienie kolumny 'LABEL' na pozycję 1 (druga kolumna)
        cols = list(labeled.columns)
        cols.insert(1, cols.pop(cols.index('LABEL')))
        labeled = labeled[cols]
    
        file_name = os.path.join(merged_dir, os.path.basename(csv_path).split('.')[0] + '_labeled.csv')
        labeled.to_csv(file_name, index=False)
    
        print(f'\nColumn {label_column_name} successfully copied from file {csv_path} to file {os.path.basename(output_path)}.')
    
    except Exception as e:
        print(f"Error occured: {e}")
    
    return labeled

