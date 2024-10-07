import pandas as pd
import os

def merger(processed_dir, csv_path):
    # Ustal ścieżkę do katalogu z plikami CSV
    merging_directory = os.path.join(os.getcwd(), processed_dir)
    
    # Sprawdź, czy katalog istnieje
    if not os.path.exists(merging_directory):
        print(f"{merging_directory} does not exist")
        return

    # Pobierz listę wszystkich plików CSV w katalogu
    csv_files = [file for file in os.listdir(merging_directory) if file.endswith(".csv")]
    
    # Posortuj listę plików
    csv_files.sort()

    # Utwórz pusty DataFrame, do którego będziesz dodawać przetworzone dane
    df_merged = pd.DataFrame()

    # Przejdź przez każdy plik CSV w katalogu
    for file in csv_files:
        file_path = os.path.join(merging_directory, file)
        
        # Wczytaj plik CSV do tymczasowego DataFrame, bez nagłówka
        df_temp = pd.read_csv(file_path, header=None)
        
        # Transponuj wczytany DataFrame
        df_temp = df_temp.transpose()
        
        # Dodaj jako pierwszą kolumnę nazwę pliku (bez rozszerzenia .csv)
        filename_without_extension = os.path.splitext(file)[0]
        df_temp.insert(0, 'filename', filename_without_extension)

        # Dodaj transponowany DataFrame z dodatkową kolumną do df_merged
        df_merged = pd.concat([df_merged, df_temp], ignore_index=True)

    # Zapisz wynikowy DataFrame do pliku CSV
    
    merged_dir = os.path.join(os.getcwd(),'merged')
    if not os.path.exists(merged_dir):
        print(f"\n{merged_dir} dir has been created")
        os.makedirs(merged_dir, exist_ok=True)
    
    file_name = os.path.basename(csv_path).split('.')[0] + '_merged.csv'
    
    output_path = os.path.join(merged_dir, file_name)
    df_merged.to_csv(output_path, index=False, header=True)
    print(f"\nMerged file saved to: {output_path}")
    
    return output_path, merged_dir


