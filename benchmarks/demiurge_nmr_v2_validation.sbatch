#!/usr/bin/env bash
#SBATCH --job-name=DEMIURGE_NMR_V2_gate
#SBATCH --partition=dgx_short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
set -Eeuo pipefail

: "${DEMIURGE_PROJECT_ROOT:?required}"
: "${DEMIURGE_VALIDATION_ROOT:?required persistent output root}"
: "${DEMIURGE_SCRATCH_ROOT:?required local scratch root}"

CONDA_ROOT=${CONDA_ROOT:-/raid/soft/miniconda}
CONDA_ENV=${CONDA_ENV:-demiurge}
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ROOT}/envs/${CONDA_ENV}"
cd -- "${DEMIURGE_PROJECT_ROOT}"

RUN_ROOT="${DEMIURGE_VALIDATION_ROOT}/job_${SLURM_JOB_ID}"
mkdir -p -- "${RUN_ROOT}" "${DEMIURGE_SCRATCH_ROOT}"
export DEMIURGE_JAVA_BUILD_DIR="${DEMIURGE_SCRATCH_ROOT}/java_build_${SLURM_JOB_ID}"
mkdir -p -m 700 -- "${DEMIURGE_JAVA_BUILD_DIR}"
python install_modules.py
REFERENCE_ROOT="${DEMIURGE_PROJECT_ROOT}/validation/frozen_nmr_v2_expected"
CORPUS="${DEMIURGE_PROJECT_ROOT}/validation/nmr_v2_parity_corpus.jsonl"
VALIDATION_INPUT=${DEMIURGE_VALIDATION_INPUT:-"${DEMIURGE_PROJECT_ROOT}/validation/nmr_v2_pipeline_smoke.csv"}

python demiurge_nmr_v2_gate.py validate-standalone --reference-root "${REFERENCE_ROOT}" --corpus "${CORPUS}" --output-root "${RUN_ROOT}/demiurge_candidate" --java-threads 2 --java-heap 4G --java-lifecycle persistent

/usr/bin/time -v -o "${RUN_ROOT}/per_batch.resources.txt" python demiurge.py run --input "${VALIDATION_INPUT}" --mode total --output-root "${RUN_ROOT}/per_batch" --temp-root "${DEMIURGE_SCRATCH_ROOT}" --prep-workers 4 --java-threads 2 --java-heap 4G --batch-size 500 --java-lifecycle per-batch --retain-scientific-artifacts --backend slurm-worker
/usr/bin/time -v -o "${RUN_ROOT}/persistent.resources.txt" python demiurge.py run --input "${VALIDATION_INPUT}" --mode total --output-root "${RUN_ROOT}/persistent" --temp-root "${DEMIURGE_SCRATCH_ROOT}" --prep-workers 4 --java-threads 2 --java-heap 4G --batch-size 500 --java-lifecycle persistent --retain-scientific-artifacts --backend slurm-worker

python demiurge_nmr_v2_gate.py compare-runs --left-output "${RUN_ROOT}/per_batch" --right-output "${RUN_ROOT}/persistent"
python demiurge_nmr_v2_gate.py compare-frozen-run --reference-root "${REFERENCE_ROOT}" --candidate-output "${RUN_ROOT}/per_batch"
python demiurge_nmr_v2_gate.py compare-frozen-run --reference-root "${REFERENCE_ROOT}" --candidate-output "${RUN_ROOT}/persistent"
python demiurge_performance_gate.py --run "per_batch=${RUN_ROOT}/per_batch" --run "persistent=${RUN_ROOT}/persistent" --output-root "${RUN_ROOT}/aggregate"
printf 'PASS standalone Demiurge NMR V2 validation: %s\n' "${RUN_ROOT}"
