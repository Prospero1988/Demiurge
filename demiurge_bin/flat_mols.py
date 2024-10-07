# flat_mols.py
import os
import subprocess

def flatten_molecules(mol_directory):
    """
    Spłaszcza struktury 3D zapisane w plikach .mol poprzez ustawienie współrzędnych Z na 0
    z wykorzystaniem Open Babel z linii poleceń.
    
    Parameters:
    mol_directory (str): Ścieżka do folderu zawierającego pliki .mol.
    """
    # Sprawdź, czy katalog istnieje
    if not os.path.exists(mol_directory):
        raise FileNotFoundError(f"The directory '{mol_directory}' does not exist.")
    
    # Pobierz listę wszystkich plików .mol w katalogu
    mol_files = [file for file in os.listdir(mol_directory) if file.endswith(".mol")]
    
    if not mol_files:
        print(f"No .mol files found in the directory '{mol_directory}'.")
        return
    
    print(f"\nFound {len(mol_files)} .mol files. Flattening structures...")

    # Iteruj przez każdy plik .mol
    counter = 0
    for mol_file in mol_files:
        mol_path = os.path.join(mol_directory, mol_file)
        
        # Przygotuj ścieżkę do pliku wyjściowego - nadpisujemy oryginał
        output_path = mol_path
        
        # Polecenie obabel, które spłaszcza strukturę przez ustawienie współrzędnych Z na 0
        command = f'obabel "{mol_path}" -O "{output_path}" --gen2D'
        
        try:
            # Uruchomienie polecenia Open Babel
            subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            counter += 1
        except subprocess.CalledProcessError as e:
            print(f"Failed to flatten {mol_file}: {e}")

    print(f"\n{counter} structures have been successfully flattened.")


