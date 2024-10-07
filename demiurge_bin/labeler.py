# labeler.py

import pandas as pd
import os

def labeler(csv_path, output_path, label_column, merged_dir):
    """
    Adds labels from a source CSV file to a merged NMR DataFrame based on a specified column and saves the result as a new CSV file.

    Parameters:
    - csv_path (str): Path to the CSV file containing labels and 'MOLECULE_NAME'.
    - output_path (str): Path to the merged CSV file containing NMR data.
    - label_column (int): Column index (1-based) in the label CSV file that contains labels to be added.
    - merged_dir (str): Directory where the labeled CSV file will be saved.

    Returns:
    - labeled (pd.DataFrame): The labeled DataFrame after adding the specified labels.
    """
    try:
        # Read the CSV file containing labels and SMILES data
        labels_csv = pd.read_csv(csv_path)
    
        # Read the merged CSV file containing NMR data
        merged_nmr = pd.read_csv(output_path)
        
        # Standardize the label_column to 0-based indexing (Python uses 0-based, user input is 1-based)
        label_column = label_column - 1
    
        # Get the name of the column at the specified index
        label_column_name = labels_csv.columns[label_column]
    
        # Merge the NMR data with the labels based on 'filename' and 'MOLECULE_NAME'
        labeled = pd.merge(merged_nmr, labels_csv[['MOLECULE_NAME', label_column_name]], 
                           left_on='filename', right_on='MOLECULE_NAME', how='left')
    
        # Assign the values from `label_column_name` to the new column `LABEL`
        labeled['LABEL'] = labeled[label_column_name]
    
        # Drop the unnecessary columns after the merge (`MOLECULE_NAME` and the original label column)
        labeled.drop(columns=['MOLECULE_NAME', label_column_name], inplace=True)
    
        # Move the 'LABEL' column to the second position in the DataFrame
        cols = list(labeled.columns)
        cols.insert(1, cols.pop(cols.index('LABEL')))
        labeled = labeled[cols]
    
        # Create the file path for saving the labeled CSV
        file_name = os.path.join(merged_dir, os.path.basename(csv_path).split('.')[0] + '_labeled.csv')
        
        # Save the labeled DataFrame to a new CSV file
        labeled.to_csv(file_name, index=False)
    
        print(f'\nColumn {label_column_name} successfully copied from file {csv_path} to file {os.path.basename(output_path)}.')
    
    except Exception as e:
        # Print error message if an exception occurs during the process
        print(f"Error occurred: {e}")
    
    return labeled


