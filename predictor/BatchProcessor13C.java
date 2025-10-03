package predictor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

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

import javax.vecmath.Point3d;

/**
 * BatchProcessor13C
 * -----------------
 * Predicts 13C NMR shifts for .mol files in parallel (8 threads).
 * - Validates 3D; if --use3d and 3D invalid, silently tries to rebuild 3D (CDK ModelBuilder3D via reflection).
 * - Suppresses CDK builder spam by muting System.out/err ONLY during rebuild (protected by a global lock).
 * - If ANY exception occurs during 3D prediction for a molecule, it silently retries that molecule in 2D
 *   and fixes the counters.
 * - Progress bar + final summary (3D native / 3D rebuilt / 2D).
 */
public class BatchProcessor13C {

    // Global sync so progress/summary prints never happen while stdout/err are muted
    private static final Object MUTE_LOCK = new Object();

    // Progress counter (for the progress bar)
    private static final AtomicInteger processedFileCount = new AtomicInteger(0);

    // Summary counters (counted once per file)
    private static final AtomicInteger count3D_native  = new AtomicInteger(0);
    private static final AtomicInteger count3D_rebuilt = new AtomicInteger(0);
    private static final AtomicInteger count2D         = new AtomicInteger(0);

    // Log "builder missing" only once (concise)
    private static final AtomicBoolean builderMissingLogged = new AtomicBoolean(false);

    // ANSI colors (optional)
    private static final String ANSI_GREEN = "\033[38;5;46m";
    private static final String ANSI_RED   = "\033[31m";
    private static final String ANSI_RESET = "\033[0m";

    private static void processMolFile(File molFile, String csvFilePath, String solvent, boolean use3d) {
        try {
            // 1) Detect V2000 vs V3000
            BufferedReader br = new BufferedReader(new FileReader(molFile));
            br.readLine(); // 1
            br.readLine(); // 2
            br.readLine(); // 3
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

            // 2) Add hydrogens
            AtomUtils.addAndPlaceHydrogens(mol);

            // 3) Aromaticity
            Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            aromaticity.apply(mol);

            // 4) Decide 3D usage (silent rebuild attempt if invalid)
            boolean use3dEffective = use3d;
            boolean rebuilt = false;

            if (use3d && hasBad3D(mol)) {
                if (regenerate3DInPlaceSilently(mol)) {
                    if (!hasBad3D(mol)) {
                        use3dEffective = true;
                        rebuilt = true;
                    } else {
                        use3dEffective = false;
                    }
                } else {
                    use3dEffective = false;
                }
            }

            // Count once per file (may be adjusted if fallback 3D->2D happens later)
            if (use3dEffective) {
                if (rebuilt) count3D_rebuilt.incrementAndGet();
                else         count3D_native.incrementAndGet();
            } else {
                count2D.incrementAndGet();
            }

            // 5) Predict and write CSV with safe retry (3D -> 2D on ANY exception)
            writePredictionsWithRetry(molFile.getName(), mol, csvFilePath, solvent, use3dEffective, rebuilt);

            processedFileCount.incrementAndGet();

        } catch (Exception e) {
            synchronized (MUTE_LOCK) {
                System.err.println(ANSI_RED + "Error while processing file "
                        + molFile.getName() + ": " + e.getMessage() + ANSI_RESET);
            }
            processedFileCount.incrementAndGet(); // keep progress moving
        }
    }

    /**
     * Try predicting in the chosen mode. If any exception occurs in 3D, redo the molecule in 2D silently
     * and fix counters, then write the CSV.
     */
    private static void writePredictionsWithRetry(String fileName,
                                                  IAtomContainer mol,
                                                  String csvFilePath,
                                                  String solvent,
                                                  boolean use3dEffectiveInitial,
                                                  boolean initialWasRebuilt3D) throws Exception {
        PredictionTool predictor = new PredictionTool();

        java.util.function.Function<Boolean, ArrayList<String>> runOnce = (use3dFlag) -> {
            ArrayList<String> lines = new ArrayList<>();
            int atomCount = mol.getAtomCount();
            for (int i = 0; i < atomCount; i++) {
                IAtom curAtom = mol.getAtom(i);
                if (curAtom.getAtomicNumber() == 6) { // 13C only
                    try {
                        // IMPORTANT: catch checked exceptions HERE and wrap → lambda stays valid.
                        float[] result = predictor.predict(mol, curAtom, use3dFlag, solvent);
                        if (result != null) {
                            lines.add(String.format(Locale.US, "%.2f%n", result[1]));
                        }
                    } catch (Exception e) {
                        throw new RuntimeException(e); // wrapped for outer retry logic
                    }
                }
            }
            return lines;
        };

        try {
            // First attempt in the initially chosen mode (3D or 2D)
            ArrayList<String> lines = runOnce.apply(use3dEffectiveInitial);
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath))) {
                for (String s : lines) writer.write(s);
            }
        } catch (RuntimeException any3dEx) {
            // Any exception during 3D -> retry 2D silently
            if (use3dEffectiveInitial) {
                // fix counters: remove the 3D count we added earlier, add 2D instead
                if (initialWasRebuilt3D) {
                    count3D_rebuilt.decrementAndGet();
                } else {
                    count3D_native.decrementAndGet();
                }
                count2D.incrementAndGet();

                ArrayList<String> lines2d = runOnce.apply(false);
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath))) {
                    for (String s : lines2d) writer.write(s);
                }
            } else {
                // already 2D; propagate
                throw any3dEx;
            }
        }
    }

    private static boolean hasBad3D(IAtomContainer mol) {
        if (mol == null) return true;
        for (IAtom a : mol.atoms()) {
            Point3d p = a.getPoint3d();
            if (p == null) return true;
            if (Double.isNaN(p.x) || Double.isNaN(p.y) || Double.isNaN(p.z)) return true;
            if (Double.isInfinite(p.x) || Double.isInfinite(p.y) || Double.isInfinite(p.z)) return true;
        }
        return false;
    }

    /**
     * Rebuild 3D in-place using CDK ModelBuilder3D via reflection, with System.out/err muted ONLY
     * inside a synchronized block protected by MUTE_LOCK (so progress/summary never print while muted).
     */
    private static boolean regenerate3DInPlaceSilently(IAtomContainer mol) {
        PrintStream originalErr;
        PrintStream originalOut;
        PrintStream devNull = new PrintStream(OutputStream.nullOutputStream());
        synchronized (MUTE_LOCK) {
            originalErr = System.err;
            originalOut = System.out;
            System.setErr(devNull);
            System.setOut(devNull);
            try {
                try {
                    Class<?> mb3dClass = Class.forName("org.openscience.cdk.modeling.builder3d.ModelBuilder3D");
                    java.lang.reflect.Method getInstance =
                            mb3dClass.getMethod("getInstance", org.openscience.cdk.interfaces.IChemObjectBuilder.class);
                    Object builder3d = getInstance.invoke(null, DefaultChemObjectBuilder.getInstance());
                    java.lang.reflect.Method generate =
                            mb3dClass.getMethod("generate3DCoordinates",
                                    org.openscience.cdk.interfaces.IAtomContainer.class, boolean.class);
                    Object newMolObj = generate.invoke(builder3d, mol, Boolean.TRUE);
                    IAtomContainer newMol = (IAtomContainer) newMolObj;
                    int n = Math.min(mol.getAtomCount(), newMol.getAtomCount());
                    for (int i = 0; i < n; i++) {
                        IAtom a = mol.getAtom(i);
                        IAtom b = newMol.getAtom(i);
                        a.setPoint3d(b.getPoint3d());
                    }
                    return true;
                } catch (ClassNotFoundException e) {
                    if (builderMissingLogged.compareAndSet(false, true)) {
                        // Print concise info AFTER unmuting (below)
                    }
                    return false;
                } catch (Throwable t) {
                    return false;
                }
            } finally {
                System.setErr(originalErr);
                System.setOut(originalOut);
                devNull.close();

                if (builderMissingLogged.compareAndSet(true, false)) {
                    synchronized (MUTE_LOCK) {
                        System.err.println("Info: CDK 3D builder (cdk-builder3d) not on classpath → skipping 3D rebuild.");
                    }
                }
            }
        }
    }

    @SuppressWarnings("unused")
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

    private static void printProgress(int current, int total) {
        int barLength = 25;
        int filledLength = (int) (barLength * ((double) current / Math.max(total, 1)));

        StringBuilder bar = new StringBuilder();
        for (int i = 0; i < filledLength; i++) bar.append('█');
        for (int i = 0; i < barLength - filledLength; i++) bar.append('-');

        int percent = (int) (100.0 * current / Math.max(total, 1));

        synchronized (MUTE_LOCK) {
            System.out.print("\rProgress: |" + bar + "| " + current + "/" + total + " (" + percent + "%)");
            System.out.flush();
            if (current >= total) System.out.println(" ");
        }
    }

    public static void main(String[] args) {
        if (args.length < 2) {
            synchronized (MUTE_LOCK) {
                System.err.println(ANSI_RED
                        + "Usage: java predictor.BatchProcessor13C <inputFolder> <outputFolder> [solvent] [no3d]"
                        + ANSI_RESET);
            }
            System.exit(1);
        }

        File inputFolder = new File(args[0]);
        File outputFolder = new File(args[1]);
        String solvent = "Unreported";
        boolean use3d = true;

        if (!inputFolder.isDirectory() || !outputFolder.isDirectory()) {
            synchronized (MUTE_LOCK) {
                System.err.println(ANSI_RED + "Input or output folder is not a directory." + ANSI_RESET);
            }
            System.exit(1);
        }

        if (args.length >= 3) solvent = args[2];
        if (args.length >= 4 && "no3d".equalsIgnoreCase(args[3])) use3d = false;

        File[] molFiles = inputFolder.listFiles((dir, name) -> name.toLowerCase().endsWith(".mol"));
        if (molFiles == null || molFiles.length == 0) {
            synchronized (MUTE_LOCK) {
                System.err.println(ANSI_RED + "No .mol files found in the input folder." + ANSI_RESET);
            }
            System.exit(1);
        }

        final File outputFolderFinal = outputFolder;
        final String solventFinal = solvent;
        final boolean use3dFinal = use3d;
        final int total = molFiles.length;

        int numThreads = 8;
        synchronized (MUTE_LOCK) {
            System.out.println("Using " + numThreads + " threads (3D set to: " + use3dFinal + ")");
        }
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (File molFile : molFiles) {
            executor.submit(() -> {
                String csvFilePath = new File(outputFolderFinal,
                        molFile.getName().replace(".mol", ".csv")).getPath();

                processMolFile(molFile, csvFilePath, solventFinal, use3dFinal);

                // Update progress bar after each file
                printProgress(processedFileCount.get(), total);
            });
        }

        executor.shutdown();
        try {
            executor.awaitTermination(7, java.util.concurrent.TimeUnit.DAYS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }

        // Close progress + summary (under lock to avoid being muted)
        printProgress(processedFileCount.get(), total);

        int n3dNative  = count3D_native.get();
        int n3dRebuilt = count3D_rebuilt.get();
        int n2d        = count2D.get();

        synchronized (MUTE_LOCK) {
            System.out.println();
            System.out.println("Summary:");
            System.out.println("  3D (native):  " + n3dNative  + " molecule(s)");
            System.out.println("  3D (rebuilt): " + n3dRebuilt + " molecule(s)");
            System.out.println("  2D:           " + n2d        + " molecule(s)");

            System.out.println(ANSI_GREEN
                    + "Total number of .mol files processed for 13C NMR prediction: "
                    + processedFileCount.get() + ANSI_RESET);
        }
    }
}
