package predictor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.util.Locale;

// Zaktualizowane importy
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader; // Nowa klasa do odczytu plików MOL
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.nmrshiftdb.PredictionTool;
import org.openscience.nmrshiftdb.util.AtomUtils;

public class BatchProcessor1H {

    private static int processedFileCount = 0;

    private static void processMolFile(File molFile, String csvFilePath, String solvent, boolean use3d, int counter) {
        try {
            // Zmodyfikowany komunikat z licznikiem
            System.out.println(counter + ". Processing file: " + molFile.getName());

            // Użycie nowego czytnika plików .mol
            MDLV2000Reader mdlreader = new MDLV2000Reader(new FileReader(molFile));
            IAtomContainer mol = (IAtomContainer) mdlreader.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));

            // Dodanie atomów wodoru
            AtomUtils.addAndPlaceHydrogens(mol);

            // Nowa metoda detekcji aromatyczności
            Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            aromaticity.apply(mol);

            // Uzyskanie predyktora
            PredictionTool predictor = new PredictionTool();

            // Przygotowanie do zapisu do pliku CSV
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath))) {
                // writer.write("mean\n"); // Możesz dodać nagłówek, jeśli potrzebujesz

                // Przetwarzanie atomów
                for (int i = 0; i < mol.getAtomCount(); i++) {
                    IAtom curAtom = mol.getAtom(i);
                    float[] result = null;

                    if (curAtom.getAtomicNumber() == 1) { // dla wodoru 1, dla wegla 6
                        result = predictor.predict(mol, curAtom, use3d, solvent);
                        if (result != null) {
                            writer.write(String.format(Locale.US, "%.2f\n", result[1]));
                        }
                    }
                }
            }

            // Zwiększ licznik przetworzonych plików
            processedFileCount++;
        } catch (Exception e) {
            System.err.println("Error processing file " + molFile.getName() + ": " + e.getMessage());
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: java BatchProcessor1H <inputFolder> <outputFolder> [solvent] [no3d]");
            System.exit(1);
        }

        File inputFolder = new File(args[0]);
        File outputFolder = new File(args[1]);
        String solvent = "Unreported";
        boolean use3d = true;

        // Sprawdzenie, czy argumenty są poprawne
        if (!inputFolder.isDirectory() || !outputFolder.isDirectory()) {
            System.err.println("Input or output folder is not a directory.");
            System.exit(1);
        }

        // Obsługa rozpuszczalnika
        if (args.length >= 3) {
            solvent = args[2];
        }

        // Obsługa opcji no3d
        if (args.length >= 4 && args[3].equalsIgnoreCase("no3d")) {
            use3d = false;
        }

        File[] molFiles = inputFolder.listFiles((dir, name) -> name.endsWith(".mol"));
        if (molFiles != null) {
            int counter = 1; // Inicjalizacja licznika plików
            for (File molFile : molFiles) {
                String csvFilePath = new File(outputFolder, molFile.getName().replace(".mol", ".csv")).getPath();
                processMolFile(molFile, csvFilePath, solvent, use3d, counter);
                counter++; // Zwiększenie licznika po każdym przetworzonym pliku
            }

            // Wyświetlenie liczby przetworzonych plików
            System.out.println("\nNumber of .mol files processed into 1H NMR predictions: " + processedFileCount);
        } else {
            System.err.println("No .mol files found in the input directory.");
        }
    }
}
