package predictor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.nmrshiftdb.PredictionTool;
import org.openscience.nmrshiftdb.util.AtomUtils;

/**
 * The BatchProcessor1H class processes a batch of .mol files in parallel (8 threads)
 * to predict 1H NMR chemical shifts. It reads molecular files, applies 3D geometry
 * and chemical property detection, and writes the results to CSV files.
 */
public class BatchProcessor1H {

    // A volatile counter to track processed .mol files (thread-safe increment).
    private static volatile int processedFileCount = 0;

    // ANSI color codes for output formatting
    private static final String ANSI_GREEN = "\033[38;5;46m";  // Green
    private static final String ANSI_RED   = "\033[31m";       // Red
    private static final String ANSI_RESET = "\033[0m";        // Reset

    /**
     * Processes a single .mol file to predict 1H NMR chemical shifts and saves the results to a CSV file.
     *
     * @param molFile     The .mol file to be processed.
     * @param csvFilePath The file path where the CSV file will be saved.
     * @param solvent     The solvent used for prediction (default "Unreported").
     * @param use3d       Whether to use 3D molecular data for the prediction.
     */
    private static void processMolFile(File molFile, String csvFilePath, String solvent, boolean use3d) {
        try {
            // 1. Determine if the file is V2000 or V3000
            BufferedReader br = new BufferedReader(new FileReader(molFile));
            br.readLine(); // Skip line 1
            br.readLine(); // Skip line 2
            br.readLine(); // Skip line 3
            String line4 = br.readLine();
            br.close();

            IAtomContainer mol;

            if (line4 != null && line4.contains("V3000")) {
                MDLV3000Reader mdlreader3000 = new MDLV3000Reader(new FileReader(molFile));
                mol = mdlreader3000.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
            } else {
                MDLV2000Reader mdlreader = new MDLV2000Reader(new FileReader(molFile));
                mol = mdlreader.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
            }

            // 2. Add hydrogens
            AtomUtils.addAndPlaceHydrogens(mol);

            // 3. Skip if molecule has down, down_inverted, or up_inverted wedge bonds
            //if (containsDownOrInvertedUpWedge(mol)) {
            //    System.err.println("Skipping molecule " + molFile.getName()
            //            + " due to down or inverted up wedge bond(s).");
            //    return;
            //}

            // 4. Apply aromaticity
            Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            aromaticity.apply(mol);

            // 5. Initialize the NMRShiftDB prediction tool
            PredictionTool predictor = new PredictionTool();

            // 6. Write predicted shifts to the CSV file
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath))) {
                int atomCount = mol.getAtomCount();
                for (int i = 0; i < atomCount; i++) {
                    IAtom curAtom = mol.getAtom(i);

                    // Predict 1H NMR shifts only if this atom is hydrogen
                    if (curAtom.getAtomicNumber() == 1) {
                        float[] result = predictor.predict(mol, curAtom, use3d, solvent);
                        if (result != null) {
                            writer.write(String.format(Locale.US, "%.2f\n", result[1]));
                        }
                    }
                }
            }

            synchronized (BatchProcessor1H.class) {
                processedFileCount++;
            }

        } catch (Exception e) {
            System.err.println(ANSI_RED + "Error while processing file "
                               + molFile.getName() + ": " + e.getMessage() + ANSI_RESET);
            e.printStackTrace();
        }
    }

    /**
     * Checks if the molecule contains a "down", "down_inverted", or "up_inverted" wedge bond.
     * Only in those cases do we skip the molecule.
     *
     * @param mol The IAtomContainer to inspect.
     * @return true if the molecule has at least one bond with DOWN, DOWN_INVERTED, or UP_INVERTED stereo.
     */
    private static boolean containsDownOrInvertedUpWedge(IAtomContainer mol) {
        if (mol == null) return false;
        for (IBond bond : mol.bonds()) {
            IBond.Stereo stereo = bond.getStereo();
            if (stereo == IBond.Stereo.DOWN
                    || stereo == IBond.Stereo.DOWN_INVERTED
                    || stereo == IBond.Stereo.UP_INVERTED) {
                return true;
            }
        }
        return false;
    }

    /**
     * Prints a dynamic progress bar to indicate the processing progress.
     *
     * @param current Number of processed files so far.
     * @param total   Total number of files to process.
     */
    private static void printProgress(int current, int total) {
        int barLength = 25;
        int filledLength = (int) (barLength * ((double) current / total));

        StringBuilder bar = new StringBuilder();
        for (int i = 0; i < filledLength; i++) {
            bar.append(ANSI_GREEN).append("█").append(ANSI_RESET);
        }
        for (int i = filledLength; i < barLength; i++) {
            bar.append("-");
        }

        int percent = (int) (100.0 * current / total);

        System.out.print("\rProgress: |" + bar + "| " + current + "/" + total + " (" + percent + "%)");
        System.out.flush();

        if (current == total) {
            System.out.println(" ");
        }
    }

    /**
     * Main method for batch processing .mol files (in parallel with 8 threads).
     *
     * Command-line arguments:
     *   args[0] - the input folder containing .mol files
     *   args[1] - the output folder for CSV files
     *   args[2] (optional) - the solvent for prediction
     *   args[3] (optional) - "no3d" to disable 3D data usage
     */
    public static void main(String[] args) {
        // Check for required arguments
        if (args.length < 2) {
            System.err.println(ANSI_RED
                    + "Usage: java BatchProcessor1H <inputFolder> <outputFolder> [solvent] [no3d]"
                    + ANSI_RESET);
            System.exit(1);
        }

        // Parse input
        File inputFolder = new File(args[0]);
        File outputFolder = new File(args[1]);
        String solvent = "Unreported"; // default
        boolean use3d = true;

        if (!inputFolder.isDirectory() || !outputFolder.isDirectory()) {
            System.err.println(ANSI_RED + "Input or output folder is not a directory." + ANSI_RESET);
            System.exit(1);
        }

        if (args.length >= 3) {
            solvent = args[2];
        }
        if (args.length >= 4 && args[3].equalsIgnoreCase("no3d")) {
            use3d = false;
        }

        // Gather all .mol files
        File[] molFiles = inputFolder.listFiles((dir, name) -> name.endsWith(".mol"));
        if (molFiles == null || molFiles.length == 0) {
            System.err.println(ANSI_RED + "No .mol files found in the input folder." + ANSI_RESET);
            return;
        }

        final int totalFiles = molFiles.length;

        // Because lambdas in older Java versions require effectively final variables:
        final File finalOutputFolder = outputFolder;
        final String finalSolvent = solvent;
        final boolean finalUse3d = use3d;

        // Create a fixed thread pool of 8
        int numThreads = 8;
        System.out.println("Using " + numThreads + " threads (3D set to: " + finalUse3d + ")");
        java.util.concurrent.ExecutorService executor
                = java.util.concurrent.Executors.newFixedThreadPool(numThreads);

        // Submit tasks for each file
        for (File molFile : molFiles) {
            executor.submit(() -> {
                String csvFilePath = new File(finalOutputFolder,
                        molFile.getName().replace(".mol", ".csv")).getPath();
                processMolFile(molFile, csvFilePath, finalSolvent, finalUse3d);
            });
        }

        // Shutdown the executor and wait for tasks to finish
        executor.shutdown();
        try {
            while (!executor.awaitTermination(500, TimeUnit.MILLISECONDS)) {
                // Update progress bar in the meantime
                printProgress(processedFileCount, totalFiles);
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // Final update of the progress bar
        printProgress(processedFileCount, totalFiles);
        System.out.println();

        System.out.println(ANSI_GREEN
                + "Total number of .mol files processed for 1H NMR prediction: "
                + processedFileCount + ANSI_RESET);
    }
}
