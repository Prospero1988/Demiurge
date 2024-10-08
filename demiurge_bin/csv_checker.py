#csv_checker.py

import pandas as pd
import csv
import os

def verify_csv(file_path):
    """
    Function to verify and modify a CSV file by handling separators, decimal points, and column structure.
    Additionally, it cleans up the first column (e.g., molecule names) by removing problematic characters.

    Parameters:
    - file_path: The path to the input CSV file to be verified.

    Returns:
    - verified_file_path: Path to the newly saved CSV file with '_verified' suffix,
                          or None if there was an error in processing the file.
    """

    try:
        # Step 1: Read a sample of the file to detect the column separator
        print(f"\nStarting verification of the file: {file_path}")
        print("\nDetecting column separator...")
        with open(file_path, 'r') as file:
            sample = file.read(2048)  # Read first 2048 characters to analyze content
        
        # Check possible delimiters by counting their occurrences in the sample
        delimiter_candidates = [',', ';', '\t']  # Common delimiters: comma, semicolon, tab
        delimiter_counts = {delim: sample.count(delim) for delim in delimiter_candidates}
        separator = max(delimiter_counts, key=delimiter_counts.get)  # Choose the most common delimiter
        print(f"\nDetected column separator based on content analysis: '{separator}'")

        # Step 2: Read the file header to determine the expected number of columns
        print("\nReading header to determine the expected number of columns...")
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter=separator)  # Read the file with detected delimiter
            headers = next(reader)  # Read header row to get column names
            expected_columns = len(headers)  # Determine the expected number of columns
            print(f"\nExpected number of columns based on header: {expected_columns}")

        # Step 3: Load CSV using the detected separator and skip rows with mismatched columns
        print("\nLoading CSV file with detected separator, skipping malformed rows...")
        # Use pandas with 'on_bad_lines' to skip malformed rows that don't match the expected column count
        df = pd.read_csv(file_path, delimiter=separator, on_bad_lines='skip')
        print(f"\nLoaded CSV file with shape after filtering rows: {df.shape}")

    except Exception as e:
        print(f"\nError reading file or detecting delimiter: {e}")
        return None

    try:
        # Step 4: Replace decimal commas with decimal dots if the separator is a semicolon
        if separator == ';':
            print("\nSeparator is semicolon. Checking for decimal commas...")
            
            # Replace commas used as decimal points only in numeric cells (not entire strings)
            # Check if cells are numeric and contain commas and do not contain dots
            def is_comma_decimal(cell):
                try:
                    # Check if the cell is a number that uses comma as decimal point
                    # Convert cell to string to handle cases where the cell is a number type
                    cell_str = str(cell)
                    # Return True if cell contains comma, does not contain dot, and is numeric
                    return ',' in cell_str and '.' not in cell_str and float(cell_str.replace(',', '.'))
                except ValueError:
                    return False

            # Apply function to detect numeric cells with commas
            if df.applymap(is_comma_decimal).any().any():
                print("\nDetected commas used as decimal points. Replacing with dots...")
                df = df.applymap(lambda x: x.replace(',', '.') if isinstance(x, str) and is_comma_decimal(x) else x)
                print("\nReplaced decimal commas with dots.")
            else:
                print("\nNo decimal commas detected.")

        # Step 5: Count the number of columns and handle excess columns
        column_count = len(df.columns)
        print(f"\nNumber of columns in the CSV file: {column_count}")
        # Remove excess columns if there are more than three
        if column_count > 3:
            print(f"\nMore than three columns detected. Removing excess columns...")
            df = df.iloc[:, :3]  # Keep only the first three columns
            print(f"\nReduced to {df.shape[1]} columns.")

        # Step 6: Change separator to comma if it was semicolon
        if separator == ';':
            print("\nChanging column separator from semicolon to comma...")
            verified_file_path = file_path.replace('.csv', '_verified.csv')  # Modify file name to include '_verified'
            df.to_csv(verified_file_path, index=False, sep=',')  # Save CSV with comma as the new separator
            print(f"\nCSV file saved with comma separator at: {verified_file_path}")
        else:
            verified_file_path = file_path.replace('.csv', '_verified.csv')
            df.to_csv(verified_file_path, index=False, sep=separator)  # Save CSV with the original separator
            print(f"\nCSV file saved with original separator at: {verified_file_path}")

        # Step 7: Clean up the first column (e.g., molecule names)
        print("\nCleaning up the first column (molecule names)...")

        def clean_molecule_name(name):
            # Remove white spaces and unwanted characters
            if isinstance(name, str):
                # Remove leading and trailing whitespace
                name = name.strip()
                # Replace any whitespace inside the name (tabs, spaces)
                name = name.replace(" ", "").replace("\t", "")
                # Replace problematic characters with underscores
                for char in ['*', '&', '^', '%', '$', '@', '!', '~', '#', '(', ')', '[', ']', '{', '}', '?', '/', '\\']:
                    name = name.replace(char, "_")
                return name
            return name

        # Apply the cleaning function to the first column
        df.iloc[:, 0] = df.iloc[:, 0].apply(clean_molecule_name)
        print("\nCleaned up molecule names in the first column.")

        # Save the cleaned file again
        df.to_csv(verified_file_path, index=False, sep=',')
        print(f"\nCSV file saved after cleaning the first column at: {verified_file_path}")

    except Exception as e:
        print(f"\nError processing CSV file: {e}")
        return None

    # Step 8: Return the path to the newly saved, verified file
    print("\nVerification and modifications completed successfully.")
    return verified_file_path
