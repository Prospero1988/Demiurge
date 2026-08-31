package predictor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Base64;
import java.util.Locale;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

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

public class BatchProcessor13C {

    private static final String PROTOCOL_PREFIX = "SPECTRAPRINTS_C3";
    private static final String PROTOCOL_VERSION = "1";
    private static PrintStream protocolOutput = System.out;

    private static final Object MUTE_LOCK = new Object();

    private static final AtomicInteger processedFileCount = new AtomicInteger(0);
    private static final AtomicInteger count3D_native = new AtomicInteger(0);
    private static final AtomicInteger count3D_rebuilt = new AtomicInteger(0);
    private static final AtomicInteger count2D = new AtomicInteger(0);

    private static final AtomicBoolean builderMissingLogged = new AtomicBoolean(false);

    private static final String ANSI_GREEN = "\033[38;5;46m";
    private static final String ANSI_RED = "\033[31m";
    private static final String ANSI_RESET = "\033[0m";

    private static final String ORIGINAL_ATOM_INDEX = "ORIGINAL_ATOM_INDEX";
    private static final String UNIFIED_PROFILE_ENV = "SPECTRAPRINTS_UNIFIED_PROFILE";
    private static final String UNIFIED_BRANCH_FILE = ".spectraprints_unified_profile_branches.jsonl";
    private static final boolean UNIFIED_PROFILE_ENABLED = environmentFlag(UNIFIED_PROFILE_ENV);

    private static Diagnostics diagnostics = Diagnostics.disabled();

    private static final class Diagnostics {
        private final File outputFile;
        private final String predictorMode;
        private final long mainStartedNanos;
        private final AtomicLong predictorConstructionNanos = new AtomicLong();
        private final AtomicLong predictorConstructionCount = new AtomicLong();
        private final AtomicLong cacheResetCount = new AtomicLong();
        private final AtomicLong cacheEntriesCleared = new AtomicLong();
        private final AtomicLong directoryScanNanos = new AtomicLong();
        private final AtomicLong molReadNanos = new AtomicLong();
        private final AtomicLong indexAnnotationNanos = new AtomicLong();
        private final AtomicLong hydrogenPreparationNanos = new AtomicLong();
        private final AtomicLong aromaticityNanos = new AtomicLong();
        private final AtomicLong coordinateValidationNanos = new AtomicLong();
        private final AtomicLong moleculePreparationNanos = new AtomicLong();
        private final AtomicLong threeDNanos = new AtomicLong();
        private final AtomicLong predictionNanos = new AtomicLong();
        private final AtomicLong csvWriteNanos = new AtomicLong();
        private final AtomicLong moleculeTotalNanos = new AtomicLong();
        private final AtomicLong rebuildExceptionCount = new AtomicLong();
        private final AtomicLong direct2DRebuildFailureCount = new AtomicLong();
        private final AtomicLong prediction3DRetry2DCount = new AtomicLong();

        private Diagnostics(File outputFile, String predictorMode) {
            this.outputFile = outputFile;
            this.predictorMode = predictorMode;
            this.mainStartedNanos = System.nanoTime();
        }

        private static Diagnostics disabled() {
            return new Diagnostics(null, "thread-local");
        }

        private boolean enabled() {
            return outputFile != null;
        }

        private long start() {
            return enabled() ? System.nanoTime() : 0L;
        }

        private void addElapsed(AtomicLong destination, long startedNanos) {
            if (enabled()) {
                destination.addAndGet(System.nanoTime() - startedNanos);
            }
        }

        private double seconds(AtomicLong nanos) {
            return nanos.get() / 1_000_000_000.0;
        }

        private void write(
                int inputFiles, int threads, String lifecycleMode,
                int batchesHandled, String shutdownReason
        ) {
            if (!enabled()) {
                return;
            }
            File parent = outputFile.getAbsoluteFile().getParentFile();
            if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
                System.err.println("Warning: cannot create Java diagnostics directory " + parent);
                return;
            }
            double wallSeconds = (System.nanoTime() - mainStartedNanos) / 1_000_000_000.0;
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
                writer.write(String.format(Locale.US,
                        "{\n"
                        + "  \"schema_version\": 3,\n"
                        + "  \"nucleus\": \"13C\",\n"
                        + "  \"predictor_mode\": \"%s\",\n"
                        + "  \"java_lifecycle\": \"%s\",\n"
                        + "  \"jvm_startup_count\": 1,\n"
                        + "  \"batches_handled\": %d,\n"
                        + "  \"shutdown_reason\": \"%s\",\n"
                        + "  \"input_files\": %d,\n"
                        + "  \"processed_files\": %d,\n"
                        + "  \"threads\": %d,\n"
                        + "  \"main_wall_seconds\": %.9f,\n"
                        + "  \"predictor_construction_count\": %d,\n"
                        + "  \"predictor_construction_seconds_sum\": %.9f,\n"
                        + "  \"per_molecule_cache_reset_count\": %d,\n"
                        + "  \"cache_entries_cleared\": %d,\n"
                        + "  \"directory_scan_seconds_sum\": %.9f,\n"
                        + "  \"mol_read_seconds_sum\": %.9f,\n"
                        + "  \"index_annotation_seconds_sum\": %.9f,\n"
                        + "  \"hydrogen_preparation_seconds_sum\": %.9f,\n"
                        + "  \"aromaticity_seconds_sum\": %.9f,\n"
                        + "  \"coordinate_validation_seconds_sum\": %.9f,\n"
                        + "  \"molecule_preparation_seconds_sum\": %.9f,\n"
                        + "  \"three_d_seconds_sum\": %.9f,\n"
                        + "  \"prediction_seconds_sum\": %.9f,\n"
                        + "  \"csv_write_seconds_sum\": %.9f,\n"
                        + "  \"molecule_total_seconds_sum\": %.9f,\n"
                        + "  \"3d_native\": %d,\n"
                        + "  \"3d_rebuilt\": %d,\n"
                        + "  \"2d\": %d,\n"
                        + "  \"rebuild_exceptions\": %d,\n"
                        + "  \"direct_2d_rebuild_failures\": %d,\n"
                        + "  \"prediction_3d_to_2d_retries\": %d\n"
                        + "}\n",
                        predictorMode, lifecycleMode, batchesHandled, shutdownReason,
                        inputFiles, processedFileCount.get(), threads,
                        wallSeconds, predictorConstructionCount.get(),
                        seconds(predictorConstructionNanos), cacheResetCount.get(),
                        cacheEntriesCleared.get(), seconds(directoryScanNanos),
                        seconds(molReadNanos), seconds(indexAnnotationNanos),
                        seconds(hydrogenPreparationNanos), seconds(aromaticityNanos),
                        seconds(coordinateValidationNanos),
                        seconds(moleculePreparationNanos), seconds(threeDNanos),
                        seconds(predictionNanos), seconds(csvWriteNanos),
                        seconds(moleculeTotalNanos), count3D_native.get(),
                        count3D_rebuilt.get(), count2D.get(),
                        rebuildExceptionCount.get(), direct2DRebuildFailureCount.get(),
                        prediction3DRetry2DCount.get()));
            } catch (Exception error) {
                System.err.println("Warning: cannot write Java diagnostics " + outputFile
                        + ": " + error.getMessage());
            }
        }
    }

    private static final class BranchRecord {
        private final String internalId;
        private String initialMode = "unavailable";
        private boolean rebuildException = false;
        private boolean prediction3DRetry2D = false;
        private String finalMode = "unavailable";
        private String finalStatus = "failure";

        private BranchRecord(String fileName) {
            this.internalId = fileName.toLowerCase(Locale.ROOT).endsWith(".mol")
                    ? fileName.substring(0, fileName.length() - 4) : fileName;
        }

        private String outcome() {
            if ("failure".equals(finalStatus)) {
                return "final_failure";
            }
            if (prediction3DRetry2D) {
                return "prediction_3d_failure_to_2d_success";
            }
            if ("rebuild_failure_2d".equals(initialMode)) {
                return "rebuild_failure_to_2d_success";
            }
            if ("native_3d".equals(initialMode)) {
                return "native_3d_success";
            }
            if ("rebuilt_3d".equals(initialMode)) {
                return "rebuilt_3d_success";
            }
            return "direct_2d_success";
        }

        private String toJson() {
            return String.format(Locale.US,
                    "{\"internal_id\":\"%s\",\"nucleus\":\"13C\","
                    + "\"initial_mode\":\"%s\",\"rebuild_exception\":%s,"
                    + "\"prediction_3d_retry_2d\":%s,\"final_mode\":\"%s\","
                    + "\"final_status\":\"%s\",\"outcome\":\"%s\"}",
                    jsonEscape(internalId), initialMode,
                    Boolean.toString(rebuildException),
                    Boolean.toString(prediction3DRetry2D), finalMode, finalStatus,
                    outcome());
        }
    }

    private static final class BatchResult {
        private final int inputFiles;
        private final int processedFiles;
        private final int native3D;
        private final int rebuilt3D;
        private final int twoD;

        private BatchResult(int inputFiles, int processedFiles, int native3D, int rebuilt3D, int twoD) {
            this.inputFiles = inputFiles;
            this.processedFiles = processedFiles;
            this.native3D = native3D;
            this.rebuilt3D = rebuilt3D;
            this.twoD = twoD;
        }
    }

    private static PredictionTool constructPredictor() {
        long started = diagnostics.start();
        try {
            return new PredictionTool();
        } catch (Exception error) {
            throw new IllegalStateException(error.getMessage(), error);
        } finally {
            if (diagnostics.enabled()) {
                diagnostics.predictorConstructionCount.incrementAndGet();
                diagnostics.addElapsed(diagnostics.predictorConstructionNanos, started);
            }
        }
    }

    private static final class PredictorLifecycle {
        private final boolean reuseByThread;
        private final ThreadLocal<PredictionTool> threadPredictor;

        private PredictorLifecycle(boolean reuseByThread) {
            this.reuseByThread = reuseByThread;
            this.threadPredictor = reuseByThread
                    ? ThreadLocal.withInitial(BatchProcessor13C::constructPredictor)
                    : null;
        }

        private PredictionTool forMolecule() {
            if (!reuseByThread) {
                return constructPredictor();
            }
            PredictionTool predictor = threadPredictor.get();
            int cachedEntries = predictor.usedHoseCodes.size();
            predictor.usedHoseCodes.clear();
            if (diagnostics.enabled()) {
                diagnostics.cacheResetCount.incrementAndGet();
                diagnostics.cacheEntriesCleared.addAndGet(cachedEntries);
            }
            return predictor;
        }
    }

    private static void processMolFile(
            File molFile,
            String csvFilePath,
            String solvent,
            boolean use3d,
            PredictorLifecycle predictorLifecycle,
            ConcurrentLinkedQueue<BranchRecord> branchRecords
    ) {
        long moleculeStarted = diagnostics.start();
        BranchRecord branch = UNIFIED_PROFILE_ENABLED ? new BranchRecord(molFile.getName()) : null;
        try {
            long started = diagnostics.start();
            IAtomContainer mol = readMolFile(molFile);
            diagnostics.addElapsed(diagnostics.molReadNanos, started);

            // Zapisujemy oryginalne indeksy WSZYSTKICH atomow z wejściowego MOL-a
            long preparationStarted = diagnostics.start();
            started = diagnostics.start();
            annotateOriginalAtomIndices(mol);
            diagnostics.addElapsed(diagnostics.indexAnnotationNanos, started);

            // Przygotowanie molekuły tak jak oczekuje predyktor
            started = diagnostics.start();
            AtomUtils.addAndPlaceHydrogens(mol);
            diagnostics.addElapsed(diagnostics.hydrogenPreparationNanos, started);

            started = diagnostics.start();
            Aromaticity aromaticity = new Aromaticity(
                    ElectronDonation.cdk(),
                    Cycles.cdkAromaticSet()
            );
            aromaticity.apply(mol);
            diagnostics.addElapsed(diagnostics.aromaticityNanos, started);
            diagnostics.addElapsed(diagnostics.moleculePreparationNanos, preparationStarted);

            boolean use3dEffective = use3d;
            boolean rebuilt = false;

            boolean bad3D = false;
            if (use3d) {
                started = diagnostics.start();
                bad3D = hasBad3D(mol);
                diagnostics.addElapsed(diagnostics.coordinateValidationNanos, started);
            }
            if (use3d && bad3D) {
                started = diagnostics.start();
                if (regenerate3DInPlaceSilently(mol, branch)) {
                    long validationStarted = diagnostics.start();
                    boolean rebuiltBad3D = hasBad3D(mol);
                    diagnostics.addElapsed(diagnostics.coordinateValidationNanos, validationStarted);
                    if (!rebuiltBad3D) {
                        use3dEffective = true;
                        rebuilt = true;
                        if (branch != null) {
                            branch.initialMode = "rebuilt_3d";
                        }
                    } else {
                        use3dEffective = false;
                    }
                } else {
                    use3dEffective = false;
                }
                diagnostics.addElapsed(diagnostics.threeDNanos, started);
                if (!use3dEffective) {
                    if (diagnostics.enabled()) {
                        diagnostics.direct2DRebuildFailureCount.incrementAndGet();
                    }
                    if (branch != null) {
                        branch.initialMode = "rebuild_failure_2d";
                    }
                }
            } else if (use3d) {
                if (branch != null) {
                    branch.initialMode = "native_3d";
                }
            } else if (branch != null) {
                branch.initialMode = "configured_2d";
            }

            if (use3dEffective) {
                if (rebuilt) {
                    count3D_rebuilt.incrementAndGet();
                } else {
                    count3D_native.incrementAndGet();
                }
            } else {
                count2D.incrementAndGet();
            }

            writePredictionsWithRetry(
                    molFile.getName(),
                    mol,
                    csvFilePath,
                    solvent,
                    use3dEffective,
                    rebuilt,
                    predictorLifecycle.forMolecule(),
                    branch
            );

            if (branch != null) {
                branch.finalStatus = "success";
                branch.finalMode = branch.prediction3DRetry2D || !use3dEffective ? "2d" : "3d";
            }

            processedFileCount.incrementAndGet();

        } catch (Exception e) {
            synchronized (MUTE_LOCK) {
                System.err.println(
                        ANSI_RED + "Error while processing file "
                        + molFile.getName() + ": " + e.getMessage() + ANSI_RESET
                );
            }
            processedFileCount.incrementAndGet();
        } finally {
            if (branch != null) {
                branchRecords.add(branch);
            }
            diagnostics.addElapsed(diagnostics.moleculeTotalNanos, moleculeStarted);
        }
    }

    private static IAtomContainer readMolFile(File molFile) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(molFile));
        br.readLine();
        br.readLine();
        br.readLine();
        String line4 = br.readLine();
        br.close();

        IAtomContainer mol;
        if (line4 != null && line4.contains("V3000")) {
            MDLV3000Reader mdlreader3000 = new MDLV3000Reader(new FileReader(molFile));
            mol = mdlreader3000.read(
                    DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class)
            );
        } else {
            MDLV2000Reader mdlreader = new MDLV2000Reader(new FileReader(molFile));
            mol = mdlreader.read(
                    DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class)
            );
        }
        return mol;
    }

    private static void annotateOriginalAtomIndices(IAtomContainer mol) {
        for (int i = 0; i < mol.getAtomCount(); i++) {
            mol.getAtom(i).setProperty(ORIGINAL_ATOM_INDEX, i + 1);
        }
    }

    private static Integer getOriginalAtomIndexOrNull(IAtom atom) {
        Object value = atom.getProperty(ORIGINAL_ATOM_INDEX);
        if (value == null) {
            return null;
        }
        if (value instanceof Integer) {
            return (Integer) value;
        }
        return Integer.parseInt(value.toString());
    }

    private static Integer getParentHeavyAtomOriginalIndexOrNull(IAtomContainer mol, IAtom hydrogenAtom) {
        for (IAtom neighbor : mol.getConnectedAtomsList(hydrogenAtom)) {
            Integer atomicNumber = neighbor.getAtomicNumber();
            if (atomicNumber != null && atomicNumber != 1) {
                return getOriginalAtomIndexOrNull(neighbor);
            }
        }
        return null;
    }

    private static void writePredictionsWithRetry(
            String fileName,
            IAtomContainer mol,
            String csvFilePath,
            String solvent,
            boolean use3dEffectiveInitial,
            boolean initialWasRebuilt3D,
            PredictionTool predictor,
            BranchRecord branch
    ) throws Exception {

        java.util.function.Function<Boolean, ArrayList<String>> runOnce = (use3dFlag) -> {
            long predictionStarted = diagnostics.start();
            ArrayList<String> lines = new ArrayList<String>();
            try {
                lines.add("mol_atom_index,element,shift\n");

                for (IAtom curAtom : mol.atoms()) {
                    Integer atomicNumber = curAtom.getAtomicNumber();
                    if (atomicNumber != null && atomicNumber == 6) {
                        try {
                            float[] result = predictor.predict(mol, curAtom, use3dFlag, solvent);
                            if (result != null) {
                                Integer originalAtomIndex = getOriginalAtomIndexOrNull(curAtom);

                                // zapisujemy TYLKO te H, ktore byly obecne w wejsciowym MOL-u
                                if (originalAtomIndex == null) {
                                    continue;
                                }

                                lines.add(String.format(
                                        Locale.US,
                                        "%d,C,%.2f%n",
                                        originalAtomIndex,
                                        result[1]
                                ));
                            }
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                    }
                }
                return lines;
            } finally {
                diagnostics.addElapsed(diagnostics.predictionNanos, predictionStarted);
            }
        };

        try {
            ArrayList<String> lines = runOnce.apply(use3dEffectiveInitial);
            long csvStarted = diagnostics.start();
            BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath));
            for (String s : lines) {
                writer.write(s);
            }
            writer.close();
            diagnostics.addElapsed(diagnostics.csvWriteNanos, csvStarted);

        } catch (RuntimeException any3dEx) {
            if (use3dEffectiveInitial) {
                if (diagnostics.enabled()) {
                    diagnostics.prediction3DRetry2DCount.incrementAndGet();
                }
                if (branch != null) {
                    branch.prediction3DRetry2D = true;
                    branch.finalMode = "2d";
                }
                if (initialWasRebuilt3D) {
                    count3D_rebuilt.decrementAndGet();
                } else {
                    count3D_native.decrementAndGet();
                }
                count2D.incrementAndGet();

                ArrayList<String> lines2d = runOnce.apply(false);
                long csvStarted = diagnostics.start();
                BufferedWriter writer = new BufferedWriter(new FileWriter(csvFilePath));
                for (String s : lines2d) {
                    writer.write(s);
                }
                writer.close();
                diagnostics.addElapsed(diagnostics.csvWriteNanos, csvStarted);

            } else {
                throw any3dEx;
            }
        }
    }

    private static boolean hasBad3D(IAtomContainer mol) {
        if (mol == null) {
            return true;
        }

        for (IAtom a : mol.atoms()) {
            Point3d p = a.getPoint3d();
            if (p == null) {
                return true;
            }
            if (Double.isNaN(p.x) || Double.isNaN(p.y) || Double.isNaN(p.z)) {
                return true;
            }
            if (Double.isInfinite(p.x) || Double.isInfinite(p.y) || Double.isInfinite(p.z)) {
                return true;
            }
        }
        return false;
    }

    private static boolean regenerate3DInPlaceSilently(IAtomContainer mol, BranchRecord branch) {
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
                    Class<?> mb3dClass = Class.forName(
                            "org.openscience.cdk.modeling.builder3d.ModelBuilder3D"
                    );

                    java.lang.reflect.Method getInstance =
                            mb3dClass.getMethod(
                                    "getInstance",
                                    org.openscience.cdk.interfaces.IChemObjectBuilder.class
                            );

                    Object builder3d = getInstance.invoke(
                            null,
                            DefaultChemObjectBuilder.getInstance()
                    );

                    java.lang.reflect.Method generate =
                            mb3dClass.getMethod(
                                    "generate3DCoordinates",
                                    org.openscience.cdk.interfaces.IAtomContainer.class,
                                    boolean.class
                            );

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
                    if (diagnostics.enabled()) {
                        diagnostics.rebuildExceptionCount.incrementAndGet();
                    }
                    if (branch != null) {
                        branch.rebuildException = true;
                    }
                    if (builderMissingLogged.compareAndSet(false, true)) {
                        // info printed after unmuting
                    }
                    return false;

                } catch (Throwable t) {
                    if (diagnostics.enabled()) {
                        diagnostics.rebuildExceptionCount.incrementAndGet();
                    }
                    if (branch != null) {
                        branch.rebuildException = true;
                    }
                    return false;
                }

            } finally {
                System.setErr(originalErr);
                System.setOut(originalOut);
                devNull.close();

                if (builderMissingLogged.compareAndSet(true, false)) {
                    synchronized (MUTE_LOCK) {
                        System.err.println(
                                "Info: CDK 3D builder (cdk-builder3d) not on classpath -> skipping 3D rebuild."
                        );
                    }
                }
            }
        }
    }

    @SuppressWarnings("unused")
    private static boolean containsDownOrInvertedUpWedge(IAtomContainer mol) {
        if (mol == null) {
            return false;
        }

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

    private static BatchResult processBatch(
            File inputFolder,
            File outputFolder,
            String solvent,
            boolean use3d,
            PredictorLifecycle predictorLifecycle,
            ExecutorService executor,
            boolean showProgress
    ) throws Exception {
        long scanStarted = diagnostics.start();
        File[] molFiles = inputFolder.listFiles(
                (dir, name) -> name.toLowerCase(Locale.ROOT).endsWith(".mol")
        );
        diagnostics.addElapsed(diagnostics.directoryScanNanos, scanStarted);
        if (molFiles == null || molFiles.length == 0) {
            throw new IllegalArgumentException("No .mol files found in the input folder.");
        }
        if (!outputFolder.isDirectory()) {
            throw new IllegalArgumentException("Output folder is not a directory.");
        }

        int processedBefore = processedFileCount.get();
        int nativeBefore = count3D_native.get();
        int rebuiltBefore = count3D_rebuilt.get();
        int twoDBefore = count2D.get();
        AtomicInteger batchProcessed = new AtomicInteger(0);
        ConcurrentLinkedQueue<BranchRecord> branchRecords = UNIFIED_PROFILE_ENABLED
                ? new ConcurrentLinkedQueue<BranchRecord>() : null;
        ArrayList<Future<?>> futures = new ArrayList<Future<?>>();
        for (File molFile : molFiles) {
            futures.add(executor.submit(() -> {
                String csvFilePath = new File(
                        outputFolder,
                        molFile.getName().replace(".mol", ".csv")
                ).getPath();
                processMolFile(
                        molFile, csvFilePath, solvent, use3d,
                        predictorLifecycle, branchRecords
                );
                int current = batchProcessed.incrementAndGet();
                if (showProgress) {
                    printProgress(current, molFiles.length);
                }
            }));
        }
        for (Future<?> future : futures) {
            future.get();
        }
        if (showProgress) {
            printProgress(batchProcessed.get(), molFiles.length);
        }
        if (UNIFIED_PROFILE_ENABLED) {
            writeBranchRecords(outputFolder, branchRecords);
        }
        return new BatchResult(
                molFiles.length,
                processedFileCount.get() - processedBefore,
                count3D_native.get() - nativeBefore,
                count3D_rebuilt.get() - rebuiltBefore,
                count2D.get() - twoDBefore
        );
    }

    private static boolean environmentFlag(String name) {
        String value = System.getenv(name);
        if (value == null) {
            return false;
        }
        String normalized = value.trim().toLowerCase(Locale.ROOT);
        return "1".equals(normalized) || "true".equals(normalized)
                || "yes".equals(normalized) || "on".equals(normalized);
    }

    private static String jsonEscape(String value) {
        return value.replace("\\", "\\\\").replace("\"", "\\\"");
    }

    private static void writeBranchRecords(
            File outputFolder, ConcurrentLinkedQueue<BranchRecord> queued
    ) {
        ArrayList<BranchRecord> records = new ArrayList<BranchRecord>(queued);
        records.sort((left, right) -> left.internalId.compareTo(right.internalId));
        File output = new File(outputFolder, UNIFIED_BRANCH_FILE);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(output))) {
            for (BranchRecord record : records) {
                writer.write(record.toJson());
                writer.newLine();
            }
        } catch (Exception error) {
            System.err.println("Warning: cannot write unified 13C branch profile "
                    + output + ": " + error.getMessage());
        }
    }

    private static String decodeProtocolValue(String encoded) {
        return new String(Base64.getUrlDecoder().decode(encoded), StandardCharsets.UTF_8);
    }

    private static String encodeProtocolValue(String value) {
        return Base64.getUrlEncoder().encodeToString(value.getBytes(StandardCharsets.UTF_8));
    }

    private static void protocolMessage(String... fields) {
        protocolOutput.println(String.join("\t", fields));
        protocolOutput.flush();
    }

    private static void runPersistentService(String[] args) {
        if (args.length != 6) {
            System.err.println("Persistent service requires solvent, 3D mode, threads, diagnostics and predictor mode.");
            System.exit(2);
        }
        String solvent = args[1];
        boolean use3d = !"no3d".equalsIgnoreCase(args[2]);
        int numThreads;
        try {
            numThreads = Integer.parseInt(args[3].trim());
        } catch (NumberFormatException error) {
            System.err.println("Invalid persistent service thread count: " + args[3]);
            System.exit(2);
            return;
        }
        if (numThreads <= 0) {
            System.err.println("Persistent service thread count must be positive.");
            System.exit(2);
            return;
        }
        File diagnosticsFile = "-".equals(args[4]) ? null : new File(args[4]);
        String predictorMode = args[5].trim().toLowerCase(Locale.ROOT);
        if (!"thread-local".equals(predictorMode) && !"per-molecule".equals(predictorMode)) {
            System.err.println("Invalid persistent predictor mode: " + args[5]);
            System.exit(2);
            return;
        }

        diagnostics = new Diagnostics(diagnosticsFile, predictorMode);
        protocolOutput = System.out;
        System.setOut(System.err);
        PredictorLifecycle predictorLifecycle = new PredictorLifecycle(
                "thread-local".equals(predictorMode)
        );
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        int inputFiles = 0;
        int batchesHandled = 0;
        String shutdownReason = "eof";
        protocolMessage(PROTOCOL_PREFIX, "READY", PROTOCOL_VERSION, "13C");
        try (BufferedReader commands = new BufferedReader(
                new InputStreamReader(System.in, StandardCharsets.UTF_8))) {
            String line;
            while ((line = commands.readLine()) != null) {
                String[] fields = line.split("\\t", -1);
                if (fields.length >= 3
                        && PROTOCOL_PREFIX.equals(fields[0])
                        && "SHUTDOWN".equals(fields[1])
                        && PROTOCOL_VERSION.equals(fields[2])) {
                    shutdownReason = "normal";
                    protocolMessage(PROTOCOL_PREFIX, "BYE", PROTOCOL_VERSION);
                    break;
                }
                String batchId = fields.length >= 3 ? fields[2] : "unknown";
                if (fields.length != 5
                        || !PROTOCOL_PREFIX.equals(fields[0])
                        || !"BATCH".equals(fields[1])) {
                    protocolMessage(PROTOCOL_PREFIX, "ERROR", batchId,
                            encodeProtocolValue("Malformed BATCH command"));
                    continue;
                }
                try {
                    File inputFolder = new File(decodeProtocolValue(fields[3]));
                    File outputFolder = new File(decodeProtocolValue(fields[4]));
                    if (!inputFolder.isDirectory() || !outputFolder.isDirectory()) {
                        throw new IllegalArgumentException("Input or output folder is not a directory.");
                    }
                    BatchResult result = processBatch(
                            inputFolder, outputFolder, solvent, use3d,
                            predictorLifecycle, executor, false
                    );
                    inputFiles += result.inputFiles;
                    batchesHandled += 1;
                    diagnostics.write(inputFiles, numThreads, "persistent", batchesHandled, "running");
                    protocolMessage(
                            PROTOCOL_PREFIX, "DONE", batchId,
                            Integer.toString(result.inputFiles),
                            Integer.toString(result.processedFiles),
                            Integer.toString(result.native3D),
                            Integer.toString(result.rebuilt3D),
                            Integer.toString(result.twoD)
                    );
                } catch (Exception error) {
                    protocolMessage(PROTOCOL_PREFIX, "ERROR", batchId,
                            encodeProtocolValue(error.toString()));
                }
            }
        } catch (Exception error) {
            shutdownReason = "protocol-error";
            System.err.println("Persistent 13C service failed: " + error.getMessage());
        } finally {
            executor.shutdownNow();
            try {
                executor.awaitTermination(30, java.util.concurrent.TimeUnit.SECONDS);
            } catch (InterruptedException error) {
                Thread.currentThread().interrupt();
                shutdownReason = "interrupted";
            }
            diagnostics.write(inputFiles, numThreads, "persistent", batchesHandled, shutdownReason);
        }
    }

    private static void printProgress(int current, int total) {
        int barLength = 25;
        int filledLength = (int) (barLength * ((double) current / Math.max(total, 1)));

        StringBuilder bar = new StringBuilder();
        for (int i = 0; i < filledLength; i++) {
            bar.append('\u2588');
        }
        for (int i = 0; i < barLength - filledLength; i++) {
            bar.append('-');
        }

        int percent = (int) (100.0 * current / Math.max(total, 1));

        synchronized (MUTE_LOCK) {
            System.out.print(
                    "\rProgress: |" + bar + "| " + current + "/" + total + " (" + percent + "%)"
            );
            System.out.flush();
            if (current >= total) {
                System.out.println(" ");
            }
        }
    }

    public static void main(String[] args) {
        if (args.length > 0 && "--persistent-service".equals(args[0])) {
            runPersistentService(args);
            return;
        }
        if (args.length < 2) {
            synchronized (MUTE_LOCK) {
                System.err.println(
                        ANSI_RED
                        + "Usage: java predictor.BatchProcessor13C <inputFolder> <outputFolder> "
                        + "[solvent] [no3d|use3d] [threads] [diagnosticsFile|-] "
                        + "[thread-local|per-molecule]"
                        + ANSI_RESET
                );
            }
            System.exit(1);
        }

        File inputFolder = new File(args[0]);
        File outputFolder = new File(args[1]);
        String solvent = "Unreported";
        boolean use3d = true;
        int numThreads = 8;
        File diagnosticsFile = null;
        String predictorMode = "thread-local";

        if (!inputFolder.isDirectory() || !outputFolder.isDirectory()) {
            synchronized (MUTE_LOCK) {
                System.err.println(
                        ANSI_RED + "Input or output folder is not a directory." + ANSI_RESET
                );
            }
            System.exit(1);
        }

        if (args.length >= 3) {
            solvent = args[2];
        }
        if (args.length >= 4 && "no3d".equalsIgnoreCase(args[3])) {
            use3d = false;
        }
        if (args.length >= 5) {
            try {
                numThreads = Integer.parseInt(args[4].trim());
            } catch (NumberFormatException error) {
                System.err.println(
                        ANSI_RED + "Invalid thread count '" + args[4]
                        + "': expected a positive integer." + ANSI_RESET
                );
                System.exit(2);
                return;
            }
            if (numThreads <= 0) {
                System.err.println(
                        ANSI_RED + "Invalid thread count '" + args[4]
                        + "': value must be greater than zero." + ANSI_RESET
                );
                System.exit(2);
                return;
            }
        }
        if (args.length >= 6 && !"-".equals(args[5])) {
            diagnosticsFile = new File(args[5]);
        }
        if (args.length >= 7) {
            predictorMode = args[6].trim().toLowerCase(Locale.ROOT);
            if (!"thread-local".equals(predictorMode)
                    && !"per-molecule".equals(predictorMode)) {
                System.err.println(
                        ANSI_RED + "Invalid predictor mode '" + args[6]
                        + "': expected thread-local or per-molecule." + ANSI_RESET
                );
                System.exit(2);
                return;
            }
        }
        diagnostics = new Diagnostics(diagnosticsFile, predictorMode);

        final PredictorLifecycle predictorLifecycle = new PredictorLifecycle(
                "thread-local".equals(predictorMode)
        );

        synchronized (MUTE_LOCK) {
            System.out.println("Using " + numThreads + " threads (3D set to: " + use3d + ")");
        }

        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        BatchResult result;
        try {
            result = processBatch(
                    inputFolder, outputFolder, solvent, use3d,
                    predictorLifecycle, executor, true
            );
        } catch (Exception error) {
            System.err.println(ANSI_RED + "Java 13C batch failed: " + error + ANSI_RESET);
            executor.shutdownNow();
            System.exit(1);
            return;
        }
        executor.shutdown();
        try {
            executor.awaitTermination(7, java.util.concurrent.TimeUnit.DAYS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }

        int n3dNative = count3D_native.get();
        int n3dRebuilt = count3D_rebuilt.get();
        int n2d = count2D.get();

        synchronized (MUTE_LOCK) {
            System.out.println();
            System.out.println("Summary:");
            System.out.println("  3D (native):  " + n3dNative + " molecule(s)");
            System.out.println("  3D (rebuilt): " + n3dRebuilt + " molecule(s)");
            System.out.println("  2D:           " + n2d + " molecule(s)");

            System.out.println(
                    ANSI_GREEN
                    + "Total number of .mol files processed for 13C NMR prediction: "
                    + processedFileCount.get()
                    + ANSI_RESET
            );
        }
        diagnostics.write(result.inputFiles, numThreads, "per-batch", 1, "process-exit");
    }
}
