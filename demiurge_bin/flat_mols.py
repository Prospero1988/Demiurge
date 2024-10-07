# flat_mols.py

import os
import subprocess

def flatten_molecules(mol_directory):
    """
    Flattens 3D structures saved in .mol files by setting Z coordinates to 0
    using Open Babel from the command line.

    Parameters:
    - mol_directory (str): Path to the directory containing .mol files.
    """
    
    # Check if the specified directory exists
    if not os.path.exists(mol_directory):
        raise FileNotFoundError(f"The directory '{mol_directory}' does not exist.")
    
    # Get a list of all .mol files in the directory
    mol_files = [file for file in os.listdir(mol_directory) if file.endswith(".mol")]
    
    # If no .mol files are found, notify and exit the function
    if not mol_files:
        print(f"No .mol files found in the directory '{mol_directory}'.")
        return
    
    print(f"\nFound {len(mol_files)} .mol files. Flattening structures...")

    # Iterate through each .mol file in the directory
    counter = 0
    for mol_file in mol_files:
        mol_path = os.path.join(mol_directory, mol_file)
        
        # Prepare the output path (in this case, it overwrites the original file)
        output_path = mol_path
        
        # Command to flatten the structure by setting Z coordinates to 0 using Open Babel
        command = f'obabel "{mol_path}" -O "{output_path}" --gen2D'
        
        try:
            # Run the Open Babel command using subprocess
            subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            counter += 1
        except subprocess.CalledProcessError as e:
            print(f"Failed to flatten {mol_file}: {e}")

    print(f"\n{counter} structures have been successfully flattened.")
