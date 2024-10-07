# custom_header.py

import pandas as pd
import os

def custom_header(labeled, csv_path):
    """
    Adds custom headers to the provided DataFrame based on the number of columns and saves it as a new CSV file.

    Parameters:
    - labeled (pd.DataFrame): The DataFrame to which custom headers will be added.
    - csv_path (str): The path to the CSV file used to generate the new filename for the output file.
    
    Returns:
    - None
    """
    try:
        # Create a basic list of headers
        header_list = ["MOLECULE_NAME", "LABEL"]
        
        # Calculate the number of additional columns that need to be named
        number_of_columns = len(labeled.columns) - 2
        
        # Initialize a counter for feature naming
        counter = 1

        # Iterate and create additional headers (FEATURE_1, FEATURE_2, etc.)
        for _ in range(number_of_columns):
            customed = "FEATURE_" + str(counter)
            header_list.append(customed)
            counter += 1
        
        # Check if the number of headers matches the number of columns in the DataFrame
        if len(header_list) == len(labeled.columns):
            labeled.columns = header_list
        else:
            print(f"Error: The number of headers ({len(header_list)}) does not match the number of columns in the DataFrame ({len(labeled.columns)}).")

        # Create a directory named 'generated_inputs' to store the final CSV file
        final_dir = os.path.join(os.getcwd(), 'generated_inputs')
        if not os.path.exists(final_dir):
            os.makedirs(final_dir, exist_ok=True)

        # Create the output filename based on the original CSV filename
        file_name = os.path.basename(csv_path).split('.')[0] + '_ML_input.csv'
        file_name = os.path.join(final_dir, file_name)

        # Save the modified DataFrame with custom headers to the output CSV file
        labeled.to_csv(file_name, index=False)
        print(f"\nFinal ML INPUT file saved as: {os.path.basename(file_name)} in {final_dir}\n")

    except Exception as e:
        # Print the error message if an exception occurs
        print(f"Error occurred: {e}")
