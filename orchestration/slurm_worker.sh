#!/usr/bin/env bash
set -Eeuo pipefail

: "${DEMIURGE_PROJECT_ROOT:?DEMIURGE_PROJECT_ROOT is required}"
: "${DEMIURGE_MANIFEST:?DEMIURGE_MANIFEST is required}"
: "${DEMIURGE_ATTEMPT:?DEMIURGE_ATTEMPT is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"

cd -- "${DEMIURGE_PROJECT_ROOT}"

mapfile -d '' -t FIELDS < <(
    python demiurge_supervisor.py row \
        --manifest "${DEMIURGE_MANIFEST}" \
        --task-index "${SLURM_ARRAY_TASK_ID}" \
        --attempt "${DEMIURGE_ATTEMPT}"
)
if ((${#FIELDS[@]} != 17)); then
    echo "Invalid manifest row field count: ${#FIELDS[@]}" >&2
    exit 2
fi

INPUT=${FIELDS[0]}
OUTPUT=${FIELDS[1]}
SCRATCH_ROOT=${FIELDS[2]}
MODE=${FIELDS[3]}
LABEL_COLUMN=${FIELDS[4]}
BATCH_SIZE=${FIELDS[5]}
PREP_WORKERS=${FIELDS[6]}
JAVA_THREADS=${FIELDS[7]}
JAVA_HEAP=${FIELDS[8]}
JAVA_LIFECYCLE=${FIELDS[9]}
MAX_ATTEMPTS=${FIELDS[10]}
CONDA_ROOT=${FIELDS[11]}
CONDA_ENV=${FIELDS[12]}
STAGING_ENABLED=${FIELDS[13]}
CAMPAIGN=${FIELDS[14]}
RETAIN_ARTIFACTS=${FIELDS[15]}
PROJECT_ROOT_FROM_MANIFEST=${FIELDS[16]}

if [[ "${PROJECT_ROOT_FROM_MANIFEST}" != "${DEMIURGE_PROJECT_ROOT}" ]]; then
    echo "Compute project root differs from frozen manifest" >&2
    exit 2
fi

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ROOT}/envs/${CONDA_ENV}"

DECISION=$(python demiurge_supervisor.py decision --manifest "${DEMIURGE_MANIFEST}" --task-index "${SLURM_ARRAY_TASK_ID}" --attempt "${DEMIURGE_ATTEMPT}")
if [[ "${DECISION}" != "RUN" ]]; then
    echo "decision=${DECISION}"
    exit 0
fi

STAGING_DIRECTORY=""
STAGING_TOKEN=""
RUNTIME_INPUT="${INPUT}"
ACTIVE_PID=""

cleanup() {
    local status=$?
    trap - EXIT INT TERM
    if [[ -n "${ACTIVE_PID}" ]]; then
        wait "${ACTIVE_PID}" 2>/dev/null || true
    fi
    if [[ -n "${STAGING_DIRECTORY}" ]]; then
        python -m orchestration.staging cleanup \
            --scratch-root "${SCRATCH_ROOT}" \
            --directory "${STAGING_DIRECTORY}" \
            --owner-token "${STAGING_TOKEN}" || {
                if ((status == 0)); then status=90; fi
            }
    fi
    python demiurge_supervisor.py record \
        --manifest "${DEMIURGE_MANIFEST}" \
        --task-index "${SLURM_ARRAY_TASK_ID}" \
        --attempt "${DEMIURGE_ATTEMPT}" \
        --exit-code "${status}" || true
    exit "${status}"
}

forward_signal() {
    local signal_name=$1
    if [[ -n "${ACTIVE_PID}" ]]; then
        kill -s "${signal_name}" "${ACTIVE_PID}" 2>/dev/null || true
    fi
}

trap cleanup EXIT
trap 'forward_signal TERM' TERM
trap 'forward_signal INT' INT

if [[ "${STAGING_ENABLED}" == "1" ]]; then
    mapfile -t STAGE_FIELDS < <(
        python -m orchestration.staging stage \
            --source "${INPUT}" \
            --scratch-root "${SCRATCH_ROOT}" \
            --campaign "${CAMPAIGN}" \
            --job-id "${SLURM_JOB_ID}" \
            --task-index "${SLURM_ARRAY_TASK_ID}" \
            --attempt "${DEMIURGE_ATTEMPT}"
    )
    STAGING_DIRECTORY=${STAGE_FIELDS[0]}
    RUNTIME_INPUT=${STAGE_FIELDS[1]}
    STAGING_TOKEN=${STAGE_FIELDS[2]}
    export DEMIURGE_JAVA_BUILD_DIR="${STAGING_DIRECTORY}/java_build"
    mkdir -p -m 700 -- "${DEMIURGE_JAVA_BUILD_DIR}"
else
    STAGING_DIRECTORY="${SCRATCH_ROOT}"
fi

mkdir -p -- "${OUTPUT}/diagnostics"
COMMAND=(python demiurge.py)
if [[ -f "${OUTPUT}/checkpoint.json" ]]; then
    COMMAND+=(resume --output-root "${OUTPUT}" --input "${RUNTIME_INPUT}")
else
    COMMAND+=(run --input "${RUNTIME_INPUT}" --canonical-input "${INPUT}" --mode "${MODE}" --output-root "${OUTPUT}" --label-column "${LABEL_COLUMN}" --max-attempts "${MAX_ATTEMPTS}")
    if [[ "${RETAIN_ARTIFACTS}" == "1" ]]; then COMMAND+=(--retain-scientific-artifacts); fi
fi
COMMAND+=(--temp-root "${STAGING_DIRECTORY}" --prep-workers "${PREP_WORKERS}" --java-threads "${JAVA_THREADS}" --java-heap "${JAVA_HEAP}" --batch-size "${BATCH_SIZE}" --java-lifecycle "${JAVA_LIFECYCLE}" --backend slurm-worker)

/usr/bin/time -v -o "${OUTPUT}/diagnostics/worker_resources_attempt_${DEMIURGE_ATTEMPT}.txt" "${COMMAND[@]}" &
ACTIVE_PID=$!
wait "${ACTIVE_PID}"
ACTIVE_PID=""
