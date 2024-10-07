# custom_header.py

import pandas as pd
import os


def custom_header(labeled, csv_path):
    try:
        # Stworzenie podstawowej listy nagłówków
        header_list = ["MOLECULE_NAME", "LABEL"]
        
        # Obliczenie liczby dodatkowych kolumn, które mają zostać nazwane
        number_of_columns = len(labeled.columns) - 2
        
        # Ustawienie licznika na początkową wartość
        counter = 1

        # Iteracja i tworzenie dodatkowych nagłówków (FEATURE_1, FEATURE_2, itd.)
        for _ in range(number_of_columns):
            customed = "FEATURE_" + str(counter)
            header_list.append(customed)
            counter += 1
        
        # Sprawdzenie, czy liczba nagłówków zgadza się z liczbą kolumn w DataFrame
        if len(header_list) == len(labeled.columns):
            labeled.columns = header_list
        else:
            print(f"Error: Liczba nagłówków ({len(header_list)}) nie zgadza się z liczbą kolumn w DataFrame ({len(labeled.columns)}).")

        final_dir = os.path.join(os.getcwd(), 'generated_inputs')
        if not os.path.exists(final_dir):
            os.makedirs(final_dir, exist_ok=True)

        file_name = os.path.basename(csv_path).split('.')[0] + '_ML_input.csv'
        file_name = os.path.join(final_dir, file_name)

        labeled.to_csv(file_name, index = False)
        print(f"\nFinal ML INPUT file saved as: {os.path.basename(file_name)} in {final_dir}\n")


    except Exception as e:
        print(f"Error occurred: {e}")