# merger.py

import pandas as pd
import os

def merger(processed_dir, csv_path):
    """
    Merges multiple CSV files in the specified directory into a single DataFrame, adds filename information, and saves the result as a new CSV file.

    Parameters:
    - processed_dir (str): Directory containing the individual CSV files to be merged.
    - csv_path (str): Path to an initial CSV file that is used to generate the name of the merged output file.

    Returns:
    - output_path (str): Path to the merged CSV file.
    - merged_dir (str): Path to the directory where the merged CSV file is saved.
    """
    
    try:
        # Determine the path to the directory with CSV files
        merging_directory = os.path.join(os.getcwd(), processed_dir)
        
        # Check if the directory exists
        if not os.path.exists(merging_directory):
            print(f"{merging_directory} does not exist.")
            return
    
        # Get a list of all CSV files in the directory
        csv_files = [file for file in os.listdir(merging_directory) if file.endswith(".csv")]
        
        # Sort the list of files
        csv_files.sort()
    
        # Create an empty DataFrame to store the merged data
        df_merged = pd.DataFrame()
    
        # Iterate through each CSV file in the directory
        for file in csv_files:
            file_path = os.path.join(merging_directory, file)
            
            # Read the CSV file into a temporary DataFrame without a header
            df_temp = pd.read_csv(file_path, header=None)
            
            # Transpose the temporary DataFrame
            df_temp = df_temp.transpose()
            
            # Add a new column at the beginning containing the filename (without the .csv extension)
            filename_without_extension = os.path.splitext(file)[0]
            df_temp.insert(0, 'filename', filename_without_extension)
    
            # Append the transposed DataFrame with the new column to the merged DataFrame
            df_merged = pd.concat([df_merged, df_temp], ignore_index=True)
    
        # Save the resulting merged DataFrame to a CSV file
        
        # Create a directory to save the merged file if it doesn't already exist
        merged_dir = os.path.join(os.getcwd(), 'merged')
        if not os.path.exists(merged_dir):
            print(f"\n{merged_dir} directory has been created.")
            os.makedirs(merged_dir, exist_ok=True)
        
        # Create a new filename for the merged file based on the original CSV filename
        file_name = os.path.basename(csv_path).split('.')[0] + '_merged.csv'
        
        # Create the full path for the merged CSV file
        output_path = os.path.join(merged_dir, file_name)
        
        # Save the merged DataFrame to the specified CSV file
        df_merged.to_csv(output_path, index=False, header=True)
        print(f"\nMerged file saved as: {os.path.basename(output_path)} in {merged_dir}")
        
    except Exception as e:
        # Print an error message if any exception occurs during the merging process
        print(f"An error occurred: {e}")
    
    return output_path, merged_dir
