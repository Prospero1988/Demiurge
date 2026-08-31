#!/usr/bin/env bash
#SBATCH --job-name=DEMIURGE_NMR_V2_gate
#SBATCH --partition=dgx_short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
set -Eeuo pipefail

: "${DEMIURGE_PROJECT_ROOT:?required}"
: "${SCREEN_PROJECT_ROOT:?required}"
: "${DEMIURGE_VALIDATION_INPUT:?required labeled CSV}"
: "${DEMIURGE_VALIDATION_ROOT:?required persistent output root}"
: "${DEMIURGE_SCRATCH_ROOT:?required local scratch root}"

CONDA_ROOT=${CONDA_ROOT:-/raid/soft/miniconda}
CONDA_ENV=${CONDA_ENV:-demiurge}
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ROOT}/envs/${CONDA_ENV}"
cd -- "${DEMIURGE_PROJECT_ROOT}"

mkdir -p -- "${DEMIURGE_VALIDATION_ROOT}" "${DEMIURGE_SCRATCH_ROOT}"
export DEMIURGE_JAVA_BUILD_DIR="${DEMIURGE_SCRATCH_ROOT}/java_build_${SLURM_JOB_ID}"
mkdir -p -m 700 -- "${DEMIURGE_JAVA_BUILD_DIR}"
python install_modules.py
python -c 'import json; from pathlib import Path; from demiurge_bin.contracts import verify_predictor_artifacts; print(json.dumps(verify_predictor_artifacts(Path.cwd()), sort_keys=True))'
CORPUS="${SCREEN_PROJECT_ROOT}/validation/nmr_v2_parity_corpus.jsonl"

python demiurge_nmr_v2_gate.py emit --implementation screen --screen-root "${SCREEN_PROJECT_ROOT}" --corpus "${CORPUS}" --output-root "${DEMIURGE_VALIDATION_ROOT}/screen_reference" --java-threads 2 --java-heap 4G --java-lifecycle persistent
python demiurge_nmr_v2_gate.py emit --implementation demiurge --corpus "${CORPUS}" --output-root "${DEMIURGE_VALIDATION_ROOT}/demiurge_candidate" --java-threads 2 --java-heap 4G --java-lifecycle persistent
python demiurge_nmr_v2_gate.py compare --screen-output "${DEMIURGE_VALIDATION_ROOT}/screen_reference" --demiurge-output "${DEMIURGE_VALIDATION_ROOT}/demiurge_candidate"

/usr/bin/time -v -o "${DEMIURGE_VALIDATION_ROOT}/per_batch.resources.txt" python demiurge.py run --input "${DEMIURGE_VALIDATION_INPUT}" --mode total --output-root "${DEMIURGE_VALIDATION_ROOT}/per_batch" --temp-root "${DEMIURGE_SCRATCH_ROOT}" --prep-workers 4 --java-threads 2 --java-heap 4G --batch-size 500 --java-lifecycle per-batch --retain-scientific-artifacts --backend slurm-worker
/usr/bin/time -v -o "${DEMIURGE_VALIDATION_ROOT}/persistent.resources.txt" python demiurge.py run --input "${DEMIURGE_VALIDATION_INPUT}" --mode total --output-root "${DEMIURGE_VALIDATION_ROOT}/persistent" --temp-root "${DEMIURGE_SCRATCH_ROOT}" --prep-workers 4 --java-threads 2 --java-heap 4G --batch-size 500 --java-lifecycle persistent --retain-scientific-artifacts --backend slurm-worker

python demiurge_nmr_v2_gate.py compare-runs --left-output "${DEMIURGE_VALIDATION_ROOT}/per_batch" --right-output "${DEMIURGE_VALIDATION_ROOT}/persistent"
python demiurge_performance_gate.py --run "per_batch=${DEMIURGE_VALIDATION_ROOT}/per_batch" --run "persistent=${DEMIURGE_VALIDATION_ROOT}/persistent" --output-root "${DEMIURGE_VALIDATION_ROOT}/aggregate"
