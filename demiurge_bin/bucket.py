# bucket.py

import os
import csv
import argparse
import numpy as np

def bucket(directory, predictor):
    """
    Function to bucket NMR spectra data based on the type of predictor (1H or 13C).
    
    Parameters:
    - directory: Directory containing input CSV files with spectra data.
    - predictor: Type of NMR predictor ('1H' or '13C').
    
    Returns:
    - processed_dir: Directory path where the bucketed spectra files are stored.
    """

    # Set the spectral window range based on the predictor type
    if predictor == "1H":
        sw_min = -1     # Minimum range for 1H NMR
        sw_max = 13     # Maximum range for 1H NMR
    elif predictor == "13C":
        sw_min = -10    # Minimum range for 13C NMR
        sw_max = 230    # Maximum range for 13C NMR

    def create_buckets():
        """
        Create 500 buckets within the range defined by predictor parameter.
        
        Returns:
        - bucket_range: Numpy array defining bucket boundaries.
        - buckets: Dictionary with 500 buckets initialized to 0.
        """
        bucket_range = np.linspace(sw_min, sw_max, 501)  # 501 points define 500 bucket intervals
        buckets = {i: 0 for i in range(500)}  # Initialize each of the 500 buckets with a count of 0
        return bucket_range, buckets
    
    def find_bucket_index(value, bucket_range):
        """
        Find the bucket index for a given spectral value.

        Parameters:
        - value: The spectral value to be assigned to a bucket.
        - bucket_range: Array defining bucket boundaries.

        Returns:
        - The index of the bucket where the value falls, or None if it is out of range.
        """
        if value < sw_min or value > sw_max:
            return None  # Value is out of the defined range
        # Iterate through the bucket range to find the appropriate index
        for i in range(len(bucket_range) - 1):
            if bucket_range[i] <= value < bucket_range[i + 1]:
                return i
        return len(bucket_range) - 2  # Value falls into the last bucket
    
    def process_file(file_path, output_dir, bucket_range):
        """
        Process a single CSV file and assign values to the corresponding buckets.

        Parameters:
        - file_path: Path to the input CSV file.
        - output_dir: Directory to store the processed file.
        - bucket_range: Array defining bucket boundaries.

        Returns:
        - error_values: List of values that couldn't be assigned to any bucket.
        """
        with open(file_path, 'r') as f:
            reader = csv.reader(f)
            buckets = {i: 0 for i in range(500)}  # Initialize buckets for this file
            error_values = []  # List to store values that cannot be assigned to buckets
            
            for row in reader:
                try:
                    # Convert the value in the row to a float and find the appropriate bucket
                    value = float(row[0])
                    bucket_index = find_bucket_index(value, bucket_range)
                    if bucket_index is not None:
                        buckets[bucket_index] += 1  # Increment the bucket count
                    else:
                        error_values.append(value)  # Out of range values
                except ValueError:
                    # If the value cannot be converted to a float, log it as an error
                    error_values.append(row[0])
            
            # Save the bucketed values to a new CSV file in the output directory
            output_file_path = os.path.join(output_dir, os.path.basename(file_path))
            with open(output_file_path, 'w') as out_f:
                for i in range(500):
                    out_f.write(f"{buckets[i]}\n")  # Write each bucket count to a new line
            
            return error_values  # Return any error values found during processing

    # Create bucket boundaries and initialize buckets
    bucket_range, _ = create_buckets()
    
    # Create a directory to store the processed files
    processed_dir = os.path.join(os.getcwd(), f'bucketed_{predictor}_spectra')
    os.makedirs(processed_dir, exist_ok=True)
    
    # Counters for success and errors
    success_count = 0
    error_files = {}  # Dictionary to track files and their errors
    
    # Process each CSV file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            file_path = os.path.join(directory, filename)
            error_values = process_file(file_path, processed_dir, bucket_range)
            
            if error_values:
                # Store errors for this file
                error_files[filename] = error_values
            else:
                # Increment success count if no errors
                success_count += 1
    
    # Print summary of the operation
    print(f"\nSuccessfully created {success_count} files as the pseudo {predictor} NMR spectra.")
    if error_files:
        print(f"Files with errors: {len(error_files)}")
        for fname, errors in error_files.items():
            print(f"{fname}: {errors}")

    return processed_dir  # Return the directory where processed files are stored
