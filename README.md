# Demiurge

Demiurge generates labeled molecular feature matrices for QSPR work. It is a standalone application: production execution, installation and validation require only this repository and its Conda environment. The production NMR path implements the validated `SPECTRAPRINTS_NMR_V2` representation historically defined with `screen_SPECTRAprints`; that repository is a reference source, not a runtime or deployment dependency. The former Demiurge ETKDG/CoordGen/OpenBabel representation is intentionally not preserved as a scientific parity target.

The project supports two first-class execution modes over one scientific core:

- local mode for hundreds to thousands of molecules on a normal Ubuntu workstation;
- SLURM mode for sharded, unattended cluster campaigns with staging, bounded retry and durable progress.

Neither backend changes the representation. Scheduler, paths, scratch location, batch size, process count, Java thread count, heap and JVM lifecycle are operational settings and are excluded from scientific identity.

## Scientific contracts

The NMR V2 path is frozen as:

```text
raw SMILES
  -> RDKit canonical isomeric SMILES
  -> explicit hydrogens
  -> Compute2DCoords
  -> RemoveStereochemistry
  -> V3000 MOL bytes
  -> Java/CDK ModelBuilder3D rebuild
  -> 3D-first NMRshiftDB2 1H and 13C prediction
  -> unnormalised count buckets
```

The active path does not call ETKDG, CoordGen or OpenBabel. Java reports whether each prediction used native 3D, rebuilt 3D or the 2D branch. Prediction JARs are SHA-256 pinned and checked before NMR work:

| Artifact | SHA-256 |
|---|---|
| `cdk-2.9.jar` | `60710218b8f9fd206e6151122e630c281462e9588e4b7a279c49c1532a8aeffe` |
| `cdk-builder3d-2.9.jar` | `2c3add480bc7363b5fe6da076f873543b47355630149927b9420126af78542ca` |
| `predictorc.jar` | `e3c3365fb3ffdccd79bb1c39c457c2486e6170f88eeaca5f36c09587950a5090` |
| `predictorh.jar` | `529e2c89279aaafcf63347460775693d0ff5120d17dae051e31b9e55f6d1e67d` |

Feature contracts and exact orders are:

| Mode | Contract ID | Exact feature order | Dimension |
|---|---|---|---:|
| `1H` | `DEMIURGE_1H_NMR_V2` | 1H buckets, `[-1, 17]`, inclusive maximum | 200 |
| `13C` | `DEMIURGE_13C_NMR_V2` | 13C buckets, `[-10, 230]`, inclusive maximum | 200 |
| `hybrid` | `DEMIURGE_HYBRID_NMR_V2_H_C` | 1H then 13C | 400 |
| `FP` | `DEMIURGE_ECFP4` | Morgan/ECFP4, radius 2, no chirality | 2048 |
| `total` | `DEMIURGE_TOTAL_NMR_V2_H_C_ECFP4` | 1H then 13C then ECFP4 | 2448 |

The migration is intentionally incompatible with models trained on the legacy Demiurge NMR representation. Such models must not consume NMR V2 features without retraining and an explicit deployment contract update. The ECFP4 component retains Demiurge's established 2048-bit contract.

## Architecture

```text
demiurge.py (local CLI) -------------------+
                                             -> demiurge_bin/pipeline.py
demiurge_supervisor.py -> SLURM worker -----+      | preparation.py
                                                    | predictor.py -> persistent Java JVMs
                                                    | bucketing.py
                                                    | atomic batch commits/run_state.py
```

`demiurge_bin/pipeline.py` is the only production feature pipeline. SLURM adds shard discovery, arrays, staging and resource policy; it does not contain molecule preparation, NMR, bucketing or fingerprint code.

Important files:

- `demiurge.py`: `run`, `resume` and read-only `status` for one input;
- `demiurge_bin/`: shared preparation, Java launcher, bucketing, composition and durable state;
- `demiurge_supervisor.py`: campaign manifest, array submission/resume and campaign status;
- `orchestration/slurm_worker.sh`: Bash worker using the same `demiurge.py` entry point;
- `orchestration/staging.py`: hash-verified, marker-owned input/scratch staging;
- `demiurge_nmr_v2_gate.py`: mandatory frozen-reference validation plus optional development cross-repository and backend/lifecycle comparison;
- `demiurge_performance_gate.py`: QC-gated performance aggregation;
- `validation/frozen_nmr_v2_expected/`: hash-pinned exact NMR V2 and total H|C|ECFP4 reference artifacts;
- `demiurge_bin/legacy_gen_mols_etkdg.py` and `demiurge-old.py`: historical reference only.

## Environment

Create or update the self-contained Python/JDK environment:

```bash
conda env create -f conda_environment.yml
# Existing environment:
conda env update -n demiurge -f conda_environment.yml --prune
conda activate demiurge
python install_modules.py
```

The environment pins OpenJDK 23.0.2 and therefore provides both `java` and `javac`; no system JDK is required. OpenBabel is deliberately absent. `python install_modules.py` reports the resolved Java/Javac paths and versions, verifies every scientific JAR hash, and proves the external build directory writable before expensive execution. Java compilation is hash-aware and stored outside the checkout. The default cache is isolated by user/job/process under `SLURM_TMPDIR` or the system temporary directory; `DEMIURGE_JAVA_BUILD_DIR` may select another owned writable directory. `DEMIURGE_JAVA` and `DEMIURGE_JAVAC` may point to explicit tools when needed. Initialization and artifact-validation errors are fatal and retain their original exception.

Input is a CSV containing `MOLECULE_NAME`, `SMILES` and the label column. `--label-column` is one-based and defaults to 3. Invalid molecules are retained as explicit failure metadata; only successful rows enter the final feature CSV.

## Local Ubuntu workflow

Run a small `total` job:

```bash
python demiurge.py run \
  --input dataset.csv \
  --mode total \
  --output-root ./results/run_001 \
  --temp-root /tmp/demiurge_run_001 \
  --prep-workers 4 \
  --java-threads 2 \
  --java-heap 4G \
  --batch-size 500 \
  --java-lifecycle persistent
```

Resume a compatible interrupted/failed run and inspect durable progress:

```bash
python demiurge.py resume --output-root ./results/run_001 --temp-root /tmp/demiurge_run_001
python demiurge.py status --output-root ./results/run_001
```

`per-batch` remains an explicit A/B/debug fallback through `--java-lifecycle per-batch`. Persistent is the production default.

## DGX/SLURM workflow

`orchestration/config.toml` contains operational defaults only. The supplied production profile is persistent JVM, Java threads 2, Java heap 4G, preparation workers 4, batch 1000, 6 CPUs and 16G per task, with at most three attempts. Adapt partition, time, environment and paths for the deployment.

```bash
python demiurge_supervisor.py submit \
  --project-root /raid/homes/$USER/Demiurge \
  --input-dir /raid/data/demiurge_shards \
  --pattern '*.csv' \
  --output-root /raid/results/demiurge \
  --scratch-root /nvme/scratch/$USER/demiurge \
  --campaign production_001

python demiurge_supervisor.py status \
  --manifest /raid/results/demiurge/production_001/campaign_manifest.json

python demiurge_supervisor.py resume \
  --manifest /raid/results/demiurge/production_001/campaign_manifest.json
```

Submission creates the stdout/stderr directory before `sbatch`. The worker is a real Bash script with `set -Eeuo pipefail`; no `--wrap`, Git checkout, NAS mount or implicit working directory is required. Inputs may be staged into a per-job scratch directory and verified by content hash. Scientific temporary MOL/raw spectra live in owned scratch; persistent batches, checkpoints, summaries, failures, diagnostics and logs remain under the output campaign. Cleanup can remove only a marker-owned child of the declared scratch root.

Three `afterany` arrays implement finite attempts. Completed/permanent tasks skip later arrays; classified transient failures may retry; unknown errors remain fail-closed. Resume refuses incompatible input content or scientific configuration and refuses while a recorded array is active. Historical job IDs absent from `squeue` are inactive; unexpected scheduler errors block resume.

## Durable outputs and recovery

Each run writes:

- `run_manifest.json`: frozen input/scientific identity and initial operational configuration;
- `checkpoint.json`: atomic state, next row, committed batches, success/failure counts and heartbeat;
- `batches/batch_*/`: atomic feature and metadata commits, optionally raw scientific artifacts;
- `generated_ML_inputs/*_ML_input.csv`: assembled final matrix;
- `failures.jsonl`: molecule-level identity, stage, type and message;
- `summary.json`: final counts, hashes, timings, throughput and observability;
- `production_progress.txt`: read-only-derived progress snapshot;
- `diagnostics/`: Java and process diagnostics.

Batch directories become visible only after their feature, metadata and commit files are complete. Resume begins at the next durable row; scratch paths and execution backend are not resume identity. A status command reads manifests/checkpoints rather than parsing stdout.

## Scientific validation

The mandatory Phase 0 deployment reference is committed under `validation/frozen_nmr_v2_expected`. Its manifest pins every expected artifact by SHA-256. These exact bytes were imported from canonical DGX job **325965** (`nmr_v2_spectraenv_20260831T130102Z`), which passed screen-reference versus Demiurge and per-batch versus persistent comparisons exactly; they were not regenerated on Windows. The originating `screen_SPECTRAprints` commit is `5c8537eeb8486b3287f0fe67c451f82f1a1a0dda`; the 12-molecule corpus SHA-256 is `fe5803b2da1356e364224e90c0cb0fc543165e4b130b38d59a4b040528b29d3e`. Normal validation verifies this inventory and never regenerates it.

Only deterministic scientific artifacts are frozen: canonical identities, preparation/prediction failure data, V3000 MOL bytes, raw indexed spectra, branch outcomes, bucket/H|C vectors, and the total H|C|ECFP4 smoke output. Checkpoints, summaries, timestamps, resource files, GC/JFR/launcher diagnostics, SLURM logs and performance aggregation are explicitly excluded. The reference is DGX/Linux byte-canonical: exact validation on a different platform may intentionally expose serializer, 2D-layout or newline differences instead of weakening the comparison. The complete old/new artifact inventory for the canonical import is recorded in `validation/nmr_v2_dgx_325965_import_manifest.json`.

The gate is zero-tolerance and fail-closed. It compares canonical identities, preparation success/failure, exact V3000 MOL bytes, indexed raw 1H/13C CSV bytes, per-molecule native-3D/rebuilt-3D/2D status, both 200-bin vectors and H|C. Full-run comparison additionally requires identical final rows, ECFP4/total composition, failures, canonical metadata and retained scientific artifacts.

Standalone validation:

```bash
python install_modules.py
python demiurge_nmr_v2_gate.py validate-standalone \
  --reference-root validation/frozen_nmr_v2_expected \
  --corpus validation/nmr_v2_parity_corpus.jsonl \
  --output-root validation/results/demiurge_candidate \
  --java-threads 2 --java-heap 4G --java-lifecycle persistent
```

The supplied DGX SBATCH additionally executes `total` in per-batch and persistent modes, compares both runs exactly, and compares each against the frozen H|C|ECFP4 result. Backend and lifecycle equivalence is tested with `compare-runs` after retaining scientific artifacts. Cross-repository emission remains available only as an optional development check when the historical checkout is present. A production recommendation is permitted only for an exact/QC-clean run.

```bash
mkdir -p /raid/homes/$USER/demiurge_validation/logs /raid/homes/$USER/demiurge_validation/runs
sbatch \
  --output=/raid/homes/$USER/demiurge_validation/logs/nmr_v2_%j.out \
  --error=/raid/homes/$USER/demiurge_validation/logs/nmr_v2_%j.err \
  --export=ALL,DEMIURGE_PROJECT_ROOT=/raid/homes/$USER/Demiurge,DEMIURGE_VALIDATION_ROOT=/raid/homes/$USER/demiurge_validation/runs,DEMIURGE_SCRATCH_ROOT=/nvme/scratch/$USER/demiurge_validation \
  /raid/homes/$USER/Demiurge/benchmarks/demiurge_nmr_v2_validation.sbatch
```

## Migration and optimization history

| Phase/candidate | Change | Validation | Status |
|---|---|---|---|
| Legacy Demiurge | ETKDG retries, CoordGen/OpenBabel fallbacks and per-batch predictor | Historical behavior; not an NMR V2 parity target | ARCHIVED |
| Phase 0 | Frozen corpus, reference provenance and exact cross-repository gate | Local official-corpus comparison against validated screen implementation: exact PASS | ACCEPTED |
| Phase 1 | Exact NMR V2 RDKit preparation and 200+200 buckets | Canonical SMILES and V3000 MOL byte parity; bucket boundary tests | ACCEPTED |
| Phase 2 | Thread-confined `PredictionTool` and per-molecule `usedHoseCodes` reset | Exact per-batch/persistent, repeated and multi-thread raw spectra | ACCEPTED |
| Phase 3 | Persistent 1H/13C JVM services with safe argv launcher and diagnostics | Exact Java integration and full-run lifecycle comparison | ACCEPTED |
| Phase 4 | Persistent preparation pool, batches, atomic commits, retry/resume | Failure/recovery and backend-equivalence tests | ACCEPTED |
| Phase 5 | Manifest-driven SLURM arrays, staging, status and bounded retries | Local orchestration tests and Bash validation; real DGX execution pending | ACCEPTED LOCALLY |
| Phase 6 | Production defaults, dependency cleanup and documentation | OpenBabel removed; hash-pinned artifacts and full local suite | ACCEPTED LOCALLY |

The corrected migration rejected two ideas: preserving legacy spectra as the parity target, and maintaining ETKDG/OpenBabel as a production/fallback mode. The known legacy 1H differences are expected evidence of the deliberate representation change, not a regression. No scientific optimization that failed NMR V2 parity was adopted.

The validated screen pipeline supplied the persistent/thread-local design evidence, but its screening throughput is not reported here as Demiurge performance. A four-molecule Windows integration smoke (two batches, 2 Java threads, 2 preparation workers) measured 8.728 s for per-batch JVM versus 6.114 s for persistent JVM, a 1.427x smoke-only speedup; all final rows and retained scientific artifacts were exact. This workload is too small for a production recommendation. The old README quoted approximately 6 minutes for 1H and 15 minutes for 13C per roughly 1000 molecules on an 8-core workstation; that was legacy code, hardware-unspecified and not a valid NMR V2 baseline. A production-like Demiurge DGX benchmark remains required before publishing representative before/after speedup.

## Current status

- **READY locally:** standalone scientific core, local CLI, checkpoint/resume, diagnostics and hash-pinned exact self-validation.
- **READY for controlled DGX validation:** the standalone SLURM gate, staging, retry/status and pinned Conda JDK are implemented without executing a cluster job from this repository task.
- **PENDING before production-scale use:** run the supplied DGX parity/lifecycle gate on a representative labeled dataset, review QC/failures and record measured resource/performance results.

## Citation and license

Leniak, A.; Pietruś, W.; Kurczab, R. *From NMR to AI: Fusing 1H and 13C Representations for Enhanced QSPR Modeling.* J. Chem. Inf. Model. 2025. [https://doi.org/10.1021/acs.jcim.5c01791](https://doi.org/10.1021/acs.jcim.5c01791).

The project is distributed under the MIT License. NMR prediction uses the NMRshiftDB2 predictor artifacts; verify their applicable terms for deployment.
