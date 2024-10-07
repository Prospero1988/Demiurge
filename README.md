# README

## NMR-Based Machine Learning Input Generator

This project provides a comprehensive pipeline for generating machine learning inputs based on feature space derived from <sup>1</sup>H or <sup>13</sup>C NMR spectra. The software reads a CSV file containing chemical compound names and their SMILES codes, processes the information to generate NMR spectra, and merges the results with a target property to create a final dataset suitable for machine learning applications.

### Features

- **Molecule Generation**: Converts SMILES codes into 3D molecular structures and saves them as `.mol` files using RDKit.
- **NMR Spectrum Prediction**: Predicts NMR spectra for each molecule using a custom Java-based [NMRshiftDB2](https://sourceforge.net/p/nmrshiftdb2/wiki/PredictorJars/) predictor.
- **Bucketization**: Converts predicted NMR spectra into a uniform matrix using a bucketing technique.
- **Data Merging**: Merges the bucketized spectra with property labels to form a consolidated dataset.
- **Label Insertion**: Adds a target property column to the merged dataset based on a specified label column.
- **Custom Headers**: Adds headers to the final dataset for easy identification and readability.
- **Optional Cleanup**: Deletes all intermediate files and folders to save space and reduce clutter.

### Requirements

Ensure the following software and libraries are installed:

1. **Python Libraries**:
   - `rdkit`
   - `pandas`
   - `numpy`

   Install the required Python packages using:

   ```bash
   pip install rdkit pandas numpy
   ```

2. **Java SDK**:
   - Java Development Kit (JDK) is required to compile and run the Java batch processor for NMR spectrum prediction. Make sure the `javac` and `java` commands are available in your system's PATH.

### Installation

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/your-username/nmr-ml-input-generator.git
cd nmr-ml-input-generator
```

### Directory Structure

The project is organized into the following directories and files:

- `demiurge_bin/`
  - `gen_mols.py`: Generates `.mol` files from SMILES strings.
  - `predictor.py`: Runs the Java batch processor for NMR spectrum prediction.
  - `bucket.py`: Converts spectra into a uniform bucketized matrix.
  - `merger.py`: Merges bucketized spectra into a single file.
  - `labeler.py`: Adds labels to the merged spectra file.
  - `custom_header.py`: Adds custom headers to the final merged dataset.
- `predictor/`
  - `BatchProcessor1H.java`: Java class for 1H NMR spectrum prediction.
  - `BatchProcessor13C.java`: Java class for 13C NMR spectrum prediction.
  - `predictorh.jar`: JAR file for 1H NMR spectrum prediction [NMRshiftDB2](https://sourceforge.net/p/nmrshiftdb2/wiki/PredictorJars/)
  - `predictorc.jar`: JAR file for 13C NMR spectrum prediction [NMRshiftDB2](https://sourceforge.net/p/nmrshiftdb2/wiki/PredictorJars/)
  - `cdk-2.9.jar`: CDK library required for spectrum prediction.

### Usage

To run the script, use the following command:

```bash
python Demiurge.py --csv_path <input_csv_file> --predictor <NMR_type> --label_column <column_number> [--clean]
```

#### Command Line Arguments

- `--csv_path`: **[Required]** Path to the input CSV file containing compound names and SMILES codes.
- `--predictor`: **[Required]** Specifies the type of NMR predictor to use (`1H` or `13C`).
- `--label_column`: **[Required]** The column index (1-based) in the input CSV file that contains the target property values.
- `--clean`: **[Optional]** If set, the script will delete all intermediate temporary files and folders after execution.

#### Example Usage

```bash
python Demiurge.py --csv_path test.csv --predictor 1H --label_column 3 --clean
```

In this example:
- The script will read the input CSV file `test.csv`.
- It will generate `.mol` files for each molecule based on its SMILES code.
- It will predict the 1H NMR spectra for each molecule.
- The spectra will be bucketized and merged into a single file.
- The target property values from column 3 in `test.csv` will be added as labels.
- All intermediate files and directories will be deleted after execution due to the `--clean` option.

### Input CSV Format

The input CSV file should have at least the following columns:
- `MOLECULE_NAME`: The name or identifier of the molecule.
- `SMILES`: The SMILES code of the molecule.
- `<Label>`: The property values to be modeled (must be specified in the `--label_column` parameter).

Example `test.csv`:

| MOLECULE_NAME | SMILES          | Label |
|---------------|-----------------|-------|
| Compound1     | CCCO            | 5.3   |
| Compound2     | CCC(=O)O        | 2.1   |
| Compound3     | CCN(CC)CC       | 7.8   |

### Script Workflow

1. **Step 1: Generate `.mol` Files**:
   - Reads SMILES codes from the input CSV file and generates corresponding `.mol` files using RDKit.

2. **Step 2: Predict NMR Spectra**:
   - Uses the Java-based `BatchProcessor1H` or `BatchProcessor13C` to predict NMR spectra for each molecule. Predictor is part of [NMRshiftDB2](https://sourceforge.net/p/nmrshiftdb2/wiki/PredictorJars/) database.

3. **Step 3: Bucketize Spectra**:
   - Converts the predicted spectra into a uniform bucketized matrix for easy analysis and machine learning input generation.

4. **Step 4: Merge Spectra and Labels**:
   - Merges the bucketized spectra with the specified label column from the input CSV file.

5. **Step 5: Add Custom Headers**:
   - Adds descriptive headers to the final merged CSV file, making it easier to interpret and use for machine learning tasks.

6. **Step 6: Cleanup (Optional)**:
   - Deletes all intermediate files and directories if the `--clean` option is specified.

### Troubleshooting

1. **Java Compilation Issues**:
   - Ensure that the `javac` and `java` commands are available and the Java SDK is installed.
   - If `javac` is not recognized, check the system's `PATH` variable and make sure it includes the path to the JDK `bin` directory.

2. **Missing Dependencies**:
   - Ensure that all required Python libraries (`rdkit`, `pandas`, and `numpy`) are installed.

3. **File Not Found Errors**:
   - Verify the paths to input files and directories. Ensure that the input CSV file and other necessary files are correctly specified.

4. **Memory or Performance Issues**:
   - If handling a large dataset, consider increasing the memory allocation for the Java runtime by adjusting the `-Xmx` parameter in the script.

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Contact

For any questions or issues, please contact [arkadiusz.leniak@gmail.com](mailto:arkadiusz.leniak@gmail.com).
