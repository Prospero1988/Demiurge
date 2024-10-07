import os
import subprocess
import argparse

def run_java_batch_processor(mol_directory, predictor):
    """
    Kompiluje i uruchamia BatchProcessor1H z wykorzystaniem Javy na wskazanym folderze mol_directory.
    
    Parameters:
    mol_directory (str): Ścieżka do katalogu wejściowego z plikami .mol.
    csv_output_folder (str): Ścieżka do katalogu wyjściowego dla plików CSV.
    """
    # Ustal ścieżkę do podkatalogu "mols" względem głównego skryptu
    csv_output_folder = os.path.join(os.getcwd(), f"predicted_spectra_{predictor}")
    
    # Sprawdź, czy katalog "csv_output_folder" istnieje, jeśli nie, to go utwórz
    if not os.path.exists(csv_output_folder):
        os.makedirs(csv_output_folder)
        print(f"\nCreated directory: {csv_output_folder}")
    
    if predictor == "1H":
        predictor_jar = ".\\predictor\\predictorh.jar"
        cdk_jar = ".\\predictor\\cdk-2.9.jar"
        batch_processor_java = ".\\predictor\\BatchProcessor1H.java"
        batch_processor_class = "predictor.BatchProcessor1H"
    elif predictor == "13C":
        predictor_jar = ".\\predictor\\predictorc.jar"
        cdk_jar = ".\\predictor\\cdk-2.9.jar"
        batch_processor_java = ".\\predictor\\BatchProcessor13C.java"
        batch_processor_class = "predictor.BatchProcessor13C"

    # Kompilacja pliku Java z -d, aby zapisać plik .class w odpowiedniej lokalizacji
    # Kompilacja pliku Java z dodatkowymi opcjami, aby wyciszyć ostrzeżenia
    compile_command = f'javac -classpath "{predictor_jar};{cdk_jar};." -d . -Xlint:-options -Xlint:deprecation -proc:none {batch_processor_java}'

    
    try:
        # Uruchomienie polecenia kompilacji
        subprocess.run(compile_command, shell=True, check=True)
        print(f"\nSuccessfully compiled {batch_processor_class}.")
        print(f"\nSpectra prediction in progress...\n")
    except subprocess.CalledProcessError as e:
        print(f"Failed to compile {batch_processor_java}: {e}")
        return

    # Dodanie ścieżki "./predictor" do classpath w run_command
    run_command = f'java -Xmx1g -classpath "{predictor_jar};{cdk_jar};.;./predictor" {batch_processor_class} "{mol_directory}" "{csv_output_folder}" "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)"'
    
    try:
        # Uruchomienie programu Java
        subprocess.run(run_command, shell=True, check=True)
        #print(f"Successfully ran {batch_processor_class} with input folder: {mol_directory} and output folder: {csv_output_folder}.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to run {batch_processor_class}: {e}")

    return csv_output_folder
