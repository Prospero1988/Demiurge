# Demiurge

![Demiurge logo](IMG/logo.png)

Demiurge is a production molecular-descriptor and machine-learning input generator. It converts labeled molecular structures into deterministic NMR-derived count vectors, ECFP4 fingerprints, or their concatenation. One scientific pipeline is shared by direct workstation runs and SLURM workers.

The current implementation is a redesign of the original Demiurge pipeline. The historical ETKDG, CoordGen, and OpenBabel preparation route is no longer active. Production uses the validated `SPECTRAPRINTS_NMR_V2` contract, deterministic RDKit preparation, CDK `ModelBuilder3D` reconstruction where required, thread-confined Java predictors, persistent JVM processes, bounded batching, atomic commits, checkpoints, retries, and hash-pinned regression fixtures.

## Scientific contract

Available modes and feature layouts are frozen:

| Mode | Feature layout | Width |
|---|---|---:|
| `1H` | indices 0–199: unnormalised 1H count buckets | 200 |
| `13C` | indices 0–199: unnormalised 13C count buckets | 200 |
| `hybrid` | 0–199: 1H; 200–399: 13C | 400 |
| `FP` | indices 0–2047: ECFP4 bits | 2048 |
| `total` | 0–199: 1H; 200–399: 13C; 400–2447: ECFP4 | 2448 |

The scientific contracts are `DEMIURGE_1H_NMR_V2`, `DEMIURGE_13C_NMR_V2`, `DEMIURGE_HYBRID_NMR_V2_H_C`, `DEMIURGE_ECFP4`, and `DEMIURGE_TOTAL_NMR_V2_H_C_ECFP4`. The declared feature dtype is `float32`; CSV output serializes the exact integer bucket counts and fingerprint bits. ECFP4 uses radius 2, 2048 bits, and no chirality flag.

NMR preparation is:

```text
raw SMILES
→ canonical isomeric SMILES
→ explicit hydrogens
→ RDKit Compute2DCoords
→ RemoveStereochemistry
→ V3000 MOL
→ Java/CDK ModelBuilder3D reconstruction
→ 3D-first 1H and 13C prediction
→ 200 + 200 count buckets
```

The NMR V2 representation is intentionally incompatible with models trained on the former legacy Demiurge spectra.

## Architecture

```text
demiurge.py / demiurge_supervisor.py
                 │
                 ▼
          input records and labels
                 │
                 ▼
       demiurge_bin/pipeline.py
          ├─ preparation.py
          ├─ predictor.py + Java/CDK
          ├─ bucketing.py
          ├─ ECFP4 composition
          └─ run_state.py / atomic batches
                 │
                 ▼
      CSV feature matrix + durable state
```

`demiurge.py` provides direct `run`, `resume`, and `status` commands. `demiurge_supervisor.py` adds only campaign discovery, SLURM arrays, staging, bounded retries, and campaign status. `orchestration/slurm_worker.sh` calls the same `demiurge.py` pipeline; orchestration contains no duplicate scientific implementation.

## Requirements and installation

The supported environment is described by `conda_environment.yml`:

- Python 3.12;
- NumPy, pandas, and RDKit;
- OpenJDK 23.0.2, including both `java` and `javac`;
- a Bash environment for SLURM execution;
- `/usr/bin/time` on workers for resource diagnostics.

Create the environment and run the production preflight:

```bash
conda env create -f conda_environment.yml
conda activate demiurge
python install_modules.py
```

The required predictor resources live in `predictor/`: `predictorh.jar`, `predictorc.jar`, `cdk-2.9.jar`, and `cdk-builder3d-2.9.jar`. Their SHA256 values are frozen in `demiurge_bin/contracts.py`; missing or changed artifacts fail closed. Java sources are compiled into an external, hash-aware cache. `DEMIURGE_JAVA_BUILD_DIR`, `DEMIURGE_JAVA`, and `DEMIURGE_JAVAC` may select explicit writable build storage or Java tools.

OpenBabel is not a production dependency.

## Input CSV

CSV input must contain exact columns `MOLECULE_NAME` and `SMILES`. The delimiter is detected automatically. `--label-column` is a one-based column position and defaults to 3; its values are copied unchanged to the output `LABEL` column and do not influence feature calculation. Blank identifiers, SMILES, or labels become explicit `INPUT_QC` failures.

`input_example.csv` is a small semicolon-delimited example.

## Standalone execution

Run the full H|C|ECFP4 representation without SLURM:

```bash
python demiurge.py run \
  --input input_example.csv \
  --mode total \
  --output-root ./results/example_total \
  --temp-root /tmp/demiurge_example_total \
  --label-column 3 \
  --batch-size 500 \
  --prep-workers 4 \
  --java-threads 2 \
  --java-heap 4G \
  --java-lifecycle persistent
```

Operational options do not enter scientific identity. `persistent` is the validated default JVM lifecycle; `--java-lifecycle per-batch` remains available for debugging and A/B validation. Use `--retain-scientific-artifacts` when durable MOL and raw NMR files are required; otherwise they remain in marker-owned temporary storage and are deleted after shutdown.

Resume a compatible interrupted run and inspect durable progress:

```bash
python demiurge.py resume \
  --output-root ./results/example_total \
  --temp-root /tmp/demiurge_example_total

python demiurge.py status --output-root ./results/example_total
```

Resume validates input content and scientific configuration exactly. Scratch paths, worker counts, batching, JVM lifecycle, and execution backend are operational settings and do not alter scientific identity.

## SLURM production execution

Production defaults in `orchestration/config.toml` are `total`, label column 3, batch size 1000, four preparation workers, two Java threads, 4G Java heap, persistent JVM, six CPUs, 16G RAM, and at most three attempts.

Submit one or more CSV shards:

```bash
/raid/soft/miniconda/envs/demiurge/bin/python /raid/homes/aleniak/Demiurge/demiurge_supervisor.py submit \
  --project-root /raid/homes/aleniak/Demiurge \
  --input-dir /raid/homes/aleniak/demiurge_inputs \
  --pattern '*.csv' \
  --output-root /raid/homes/aleniak/demiurge_runs \
  --scratch-root /nvme/scratch/aleniak/demiurge_runs/production_001 \
  --campaign production_001
```

Important `submit` options are:

- discovery: `--input-dir`, `--pattern`, `--campaign`;
- durable/runtime paths: `--output-root`, `--scratch-root`, `--project-root`;
- science selection: `--mode`, `--label-column`;
- operational tuning: `--batch-size`, `--prep-workers`, `--java-threads`, `--java-heap`, `--java-lifecycle`;
- scheduler resources: `--cpus-per-task`, `--memory`, `--partition`, `--time`, `--max-concurrent`, `--job-name-prefix`;
- deployment: `--conda-root`, `--conda-env`, `--no-staging`, `--retain-scientific-artifacts`.

`/nvme` is compute-node-local DGX storage. Do not create it from the login node. The path is recorded in the campaign manifest without requiring it to exist there; after SLURM starts, the worker creates an isolated marker-owned child, stages and hash-verifies the input, places temporary MOL/raw NMR/build data there, and removes only its owned directory during cleanup.

Submission creates durable stdout/stderr directories before calling `sbatch`. The worker is a real Bash script with `set -Eeuo pipefail`; users should paste only the supervisor command into an interactive login shell, not wrap interactive commands in shell-wide `set -e` or `exit` guards.

### Campaign status, retries, and resume

```bash
python demiurge_supervisor.py status \
  --manifest /raid/homes/aleniak/demiurge_runs/production_001/campaign_manifest.json

python demiurge_supervisor.py resume \
  --manifest /raid/homes/aleniak/demiurge_runs/production_001/campaign_manifest.json
```

The initial submission creates a dependency chain of at most three attempts. Each worker checks durable state before running; already completed tasks are skipped. Resume refuses to overlap an active recorded SLURM array, treats scheduler-expired historical job IDs as inactive, and remains fail-closed on unexpected scheduler errors. Status reads manifests and checkpoints rather than parsing logs and writes `production_progress.txt`.

## Output layout

For campaign `<output-root>/<campaign>/`:

```text
campaign_manifest.json
production_progress.txt
logs/<array>_<task>.{out,err}
attempt_history/task-*/attempt-*.json
results/<input-stem>_<input-sha256-prefix>/
  run_manifest.json
  checkpoint.json
  production_progress.txt
  summary.json
  failures.jsonl
  diagnostics/
  batches/batch_*/
    features.csv
    metadata.jsonl
    commit.json
    scientific_artifacts/        # only with --retain-scientific-artifacts
  generated_ML_inputs/<input-stem>_<mode>_ML_input.csv
```

Only successful molecules enter the final feature matrix. Each batch records all molecule outcomes in `metadata.jsonl`; failures are aggregated into `failures.jsonl`. Atomic batch directories are the resume boundary. The summary stores counts, timings, throughput, the final output path, and its SHA256.

## Validation and reproducibility

The hash-pinned exact reference under `validation/frozen_nmr_v2_expected/` originates from validated DGX job 325965 and was not regenerated on Windows. It contains only deterministic scientific artifacts; GC logs, process IDs, timestamps, resource telemetry, and launcher diagnostics are excluded from the frozen contract.

The validation gate compares canonical identity, V3000 MOL bytes, preparation outcomes, indexed raw 1H/13C CSVs, 3D branch status, NMR buckets, final H|C vectors, and total H|C|ECFP4 rows with zero tolerance. DGX job 326101 passed production dependency preflight, exact standalone NMR V2 parity, persistent/per-batch parity, and frozen total parity. The validated total fixture SHA256 is:

```text
cab58dfc7c6ad0f835f30e4fa1929aca4882c47c19bdcf3841355a5b0c5cf607
```

Run the standalone frozen checks with:

```bash
python demiurge_nmr_v2_gate.py verify-reference \
  --reference-root validation/frozen_nmr_v2_expected \
  --corpus validation/nmr_v2_parity_corpus.jsonl
```

The complete DGX validation wrapper is `benchmarks/demiurge_nmr_v2_validation.sbatch`. Exact MOL bytes can depend on validated RDKit/platform serialization; the DGX fixture is authoritative and comparisons intentionally remain byte-exact.

## Repository layout

- `demiurge.py` — direct production CLI;
- `demiurge_bin/` — active scientific core and durable state;
- `demiurge_supervisor.py` — production campaign orchestration;
- `orchestration/` — SLURM worker, staging, and operational defaults;
- `predictor/` — hash-pinned Java/CDK resources and Java sources;
- `validation/` — canonical corpus and frozen scientific fixtures;
- `benchmarks/` — controlled validation wrapper;
- `tests/` — scientific, lifecycle, failure, and orchestration regression tests.

The obsolete original execution stack is retained in Git history and on the `GPT_corrected` reference branch, not in the production source tree.

## Citation and license

Leniak, A.; Pietruś, W.; Kurczab, R. *From NMR to AI: Fusing 1H and 13C Representations for Enhanced QSPR Modeling.* J. Chem. Inf. Model. 2025. <https://doi.org/10.1021/acs.jcim.5c01791>.

See `LICENSE` for repository licensing terms.
