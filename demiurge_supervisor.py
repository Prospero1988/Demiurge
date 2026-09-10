#!/usr/bin/env python3
"""Optional SLURM campaign backend for the shared Demiurge NMR V2 pipeline."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any, Callable

from demiurge_bin.contracts import object_sha256, scientific_contract, verify_predictor_artifacts
from demiurge_bin.io import create_input_reader
from demiurge_bin.java_heap import normalize_java_heap
from demiurge_bin.run_state import atomic_write_json, atomic_write_text, input_identity, read_json, utc_now


MAX_ATTEMPTS = 3
Runner = Callable[..., subprocess.CompletedProcess[str]]


def positive_integer(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def column_selector(value: str) -> int | str:
    selected = value.strip()
    if not selected:
        raise argparse.ArgumentTypeError("column selector must not be empty")
    if selected.isdecimal():
        return positive_integer(selected)
    return selected


def load_defaults(project_root: Path) -> dict[str, Any]:
    with (project_root / "orchestration" / "config.toml").open("rb") as handle:
        return tomllib.load(handle)


def _safe_job_name(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.-")
    if not normalized:
        raise ValueError("SLURM job name is empty after sanitization")
    return normalized[:96]


def _manifest_contract(document: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in document.items() if key not in {"submissions", "contract_sha256"}}


def save_manifest(path: Path, document: dict[str, Any]) -> None:
    document["contract_sha256"] = object_sha256(_manifest_contract(document))
    atomic_write_json(path, document)


def load_manifest(path: Path) -> tuple[Path, dict[str, Any]]:
    selected = path.expanduser().resolve()
    document = read_json(selected)
    if int(document.get("manifest_schema_version", -1)) != 1:
        raise RuntimeError("Unsupported campaign manifest schema")
    expected = document.get("contract_sha256")
    if expected != object_sha256(_manifest_contract(document)):
        raise RuntimeError("Campaign manifest frozen contract hash mismatch")
    if int(document.get("max_attempts", -1)) != MAX_ATTEMPTS:
        raise RuntimeError(f"Demiurge campaigns require max_attempts={MAX_ATTEMPTS}")
    return selected, document


def prepare_manifest(args: argparse.Namespace) -> tuple[Path, dict[str, Any]]:
    project_root = args.project_root.expanduser().resolve()
    if not (project_root / "demiurge.py").is_file():
        raise RuntimeError(f"Invalid Demiurge project root: {project_root}")
    defaults = load_defaults(project_root)
    worker = defaults["worker"]
    slurm = defaults["slurm"]
    input_dir = args.input_dir.expanduser().resolve()
    inputs = sorted(path.resolve() for path in input_dir.glob(args.pattern) if path.is_file())
    if not inputs:
        raise RuntimeError(f"No input shards match {args.pattern!r} in {input_dir}")
    campaign_root = args.output_root.expanduser().resolve() / args.campaign
    campaign_root.mkdir(parents=True, exist_ok=True)
    manifest_path = campaign_root / "campaign_manifest.json"
    if manifest_path.exists():
        raise RuntimeError(f"Campaign manifest already exists: {manifest_path}")
    mode = args.mode or worker["mode"]
    selected_label = getattr(args, "label_column", None)
    label_column = selected_label if selected_label is not None else worker["label_column"]
    input_format = getattr(args, "input_format", "csv")
    input_table = getattr(args, "input_table", None)
    input_query = getattr(args, "input_query", None)
    id_column = getattr(args, "id_column", "MOLECULE_NAME")
    smiles_column = getattr(args, "smiles_column", "SMILES")
    output_format = getattr(args, "output_format", "csv")
    output_table = getattr(args, "output_table", "demiurge_features")
    metadata_table = getattr(args, "metadata_table", "demiurge_metadata")
    io_config = {
        "input_format": input_format,
        "input_table": input_table,
        "input_query": input_query,
        "id_column": id_column,
        "smiles_column": smiles_column,
        "output_format": output_format,
        "output_table": output_table,
        "metadata_table": metadata_table,
    }
    contract = scientific_contract(mode)
    tasks = []
    for index, path in enumerate(inputs):
        digest = input_identity(path)
        output = campaign_root / "results" / f"{path.stem}_{digest['sha256'][:12]}"
        description = create_input_reader(
            path,
            input_format=input_format,
            input_table=input_table,
            input_query=input_query,
            id_column=id_column,
            smiles_column=smiles_column,
            label_column=label_column,
        ).describe()
        tasks.append({
            "task_index": index,
            "input_path": str(path),
            "input_identity": digest,
            "expected_rows": description.total_rows,
            "output_root": str(output),
        })
    document = {
        "manifest_schema_version": 1,
        "campaign": args.campaign,
        "created_at": utc_now(),
        "project_root": str(project_root),
        "input_dir": str(input_dir),
        "pattern": args.pattern,
        "campaign_root": str(campaign_root),
        "scratch_root": str(args.scratch_root.expanduser().resolve()),
        "max_attempts": MAX_ATTEMPTS,
        "scientific": {
            "mode": mode,
            "label_column": label_column,
            "contract": contract,
            "contract_sha256": object_sha256(contract),
            "predictor_artifact_sha256": verify_predictor_artifacts(project_root) if mode != "FP" else {},
        },
        "io": io_config,
        "resources": {
            "batch_size": int(args.batch_size or worker["batch_size"]),
            "prep_workers": int(args.prep_workers or worker["prep_workers"]),
            "java_threads": int(args.java_threads or worker["java_threads"]),
            "java_heap": normalize_java_heap(args.java_heap or worker["java_heap"]),
            "java_lifecycle": args.java_lifecycle or worker["java_lifecycle"],
            "retain_scientific_artifacts": bool(args.retain_scientific_artifacts),
        },
        "slurm": {
            "partition": args.partition or slurm["partition"],
            "time": args.time or slurm["time"],
            "cpus_per_task": int(args.cpus_per_task or slurm["cpus_per_task"]),
            "memory": args.memory or slurm["memory"],
            "max_concurrent_jobs": int(args.max_concurrent or slurm["max_concurrent_jobs"]),
            "job_name_prefix": _safe_job_name(args.job_name_prefix or slurm["job_name_prefix"]),
            "conda_root": args.conda_root or slurm["conda_root"],
            "conda_env": args.conda_env or slurm["conda_env"],
            "staging_enabled": bool(slurm["staging_enabled"] and not args.no_staging),
        },
        "tasks": tasks,
        "submissions": [],
    }
    save_manifest(manifest_path, document)
    return manifest_path, document


def submit_array(
    manifest_path: Path,
    manifest: dict[str, Any],
    task_indices: list[int],
    attempt: int,
    dependency: str | None,
    *,
    runner: Runner = subprocess.run,
    dry_run: bool = False,
) -> str:
    if not task_indices:
        raise ValueError("Cannot submit an empty task array")
    campaign_root = Path(manifest["campaign_root"])
    logs = campaign_root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    slurm = manifest["slurm"]
    indices = ",".join(str(index) for index in sorted(set(task_indices)))
    array = f"{indices}%{slurm['max_concurrent_jobs']}"
    export = (
        f"ALL,DEMIURGE_PROJECT_ROOT={manifest['project_root']},"
        f"DEMIURGE_MANIFEST={manifest_path},DEMIURGE_ATTEMPT={attempt}"
    )
    command = [
        "sbatch", "--parsable",
        "--job-name", f"{slurm['job_name_prefix']}-{manifest['campaign']}-a{attempt:02d}",
        "--partition", str(slurm["partition"]),
        "--time", str(slurm["time"]),
        "--cpus-per-task", str(slurm["cpus_per_task"]),
        "--mem", str(slurm["memory"]),
        "--array", array,
        "--output", str(logs / "%A_%a.out"),
        "--error", str(logs / "%A_%a.err"),
        "--export", export,
    ]
    if dependency:
        command.extend(["--dependency", f"afterany:{dependency}"])
    command.append(str(Path(manifest["project_root"]) / "orchestration" / "slurm_worker.sh"))
    if dry_run:
        job_id = f"DRYRUN{attempt}"
    else:
        completed = runner(command, check=False, text=True, capture_output=True)
        if completed.returncode != 0:
            raise RuntimeError(f"sbatch failed ({completed.returncode}): {completed.stderr.strip()}")
        job_id = completed.stdout.strip().split(";", 1)[0]
        if not job_id:
            raise RuntimeError("sbatch returned no job ID")
    manifest["submissions"].append({
        "attempt": attempt,
        "job_id": job_id,
        "task_indices": sorted(set(task_indices)),
        "dependency": dependency,
        "submitted_at": utc_now(),
        "command": command,
    })
    save_manifest(manifest_path, manifest)
    return job_id


def submit_chain(
    manifest_path: Path,
    manifest: dict[str, Any],
    task_indices: list[int],
    *,
    start_attempt: int = 1,
    runner: Runner = subprocess.run,
    dry_run: bool = False,
) -> list[str]:
    dependency = None
    jobs = []
    for attempt in range(start_attempt, MAX_ATTEMPTS + 1):
        job = submit_array(manifest_path, manifest, task_indices, attempt, dependency, runner=runner, dry_run=dry_run)
        jobs.append(job)
        dependency = job
    return jobs


def _attempt_path(manifest_path: Path, task_index: int, attempt: int) -> Path:
    return manifest_path.parent / "attempt_history" / f"task-{task_index:04d}" / f"attempt-{attempt:02d}.json"


def classify_failure(message: str) -> tuple[str, bool]:
    lowered = message.lower()
    permanent_tokens = (
        "scientific configuration", "contract", "hash mismatch", "input content identity",
        "checkpoint schema", "missing columns", "label_column", "required scientific predictor",
    )
    transient_tokens = (
        "timeout", "timed out", "preempt", "node failure", "broken pipe", "interrupted",
        "i/o", "input/output error", "temporarily unavailable", "java batch subprocess failed",
    )
    if any(token in lowered for token in permanent_tokens):
        return "PERMANENT_CONFIGURATION", False
    if any(token in lowered for token in transient_tokens):
        return "TRANSIENT_OPERATIONAL", True
    return "UNCLASSIFIED_FAIL_CLOSED", False


def retry_decision(manifest_path: Path, manifest: dict[str, Any], task_index: int, attempt: int) -> str:
    if attempt < 1 or attempt > MAX_ATTEMPTS:
        return "SKIP_EXHAUSTED"
    task = manifest["tasks"][task_index]
    checkpoint_path = Path(task["output_root"]) / "checkpoint.json"
    if checkpoint_path.is_file():
        checkpoint = read_json(checkpoint_path)
        if checkpoint.get("status") == "DONE":
            return "SKIP_DONE"
    for prior in range(1, attempt):
        path = _attempt_path(manifest_path, task_index, prior)
        if path.is_file() and read_json(path).get("retryable") is False:
            return "SKIP_PERMANENT"
    return "RUN"


def record_attempt(manifest_path: Path, manifest: dict[str, Any], task_index: int, attempt: int, exit_code: int) -> dict[str, Any]:
    task = manifest["tasks"][task_index]
    checkpoint_path = Path(task["output_root"]) / "checkpoint.json"
    checkpoint = read_json(checkpoint_path) if checkpoint_path.is_file() else {}
    status = str(checkpoint.get("status", "MISSING"))
    message = str(checkpoint.get("failure_message") or f"worker exit code {exit_code}, checkpoint={status}")
    if exit_code == 0 and status == "DONE":
        classification, retryable, final_status = "SUCCESS", False, "DONE"
    else:
        classification, retryable = classify_failure(message)
        final_status = "FAILED_EXHAUSTED" if attempt >= MAX_ATTEMPTS else "FAILED"
        retryable = retryable and attempt < MAX_ATTEMPTS
    record = {
        "schema_version": 1,
        "task_index": task_index,
        "attempt": attempt,
        "exit_code": exit_code,
        "status": final_status,
        "checkpoint_status": status,
        "classification": classification,
        "retryable": retryable,
        "failure_message": None if final_status == "DONE" else message,
        "recorded_at": utc_now(),
    }
    atomic_write_json(_attempt_path(manifest_path, task_index, attempt), record)
    return record


def ensure_previous_arrays_inactive(manifest: dict[str, Any], runner: Runner = subprocess.run) -> None:
    active = []
    for submission in manifest.get("submissions", []):
        job_id = str(submission["job_id"])
        if job_id.startswith("DRYRUN"):
            continue
        result = runner(["squeue", "-h", "-j", job_id, "-o", "%i|%T"], check=False, text=True, capture_output=True)
        if result.returncode != 0:
            combined = f"{result.stdout}\n{result.stderr}".lower()
            if "invalid job id specified" in combined:
                continue
            raise RuntimeError(f"Could not verify previous SLURM array {job_id}: {result.stderr.strip()}")
        if result.stdout.strip():
            active.append(result.stdout.strip())
    if active:
        raise RuntimeError("Previous production array is still active:\n" + "\n".join(active))


def campaign_status(manifest_path: Path, manifest: dict[str, Any]) -> tuple[str, Path]:
    totals = {
        "expected": sum(int(task["expected_rows"]) for task in manifest["tasks"]),
        "successful": 0,
        "failed_molecules": 0,
        "completed_shards": 0,
        "running_shards": 0,
        "pending_shards": 0,
        "failed_shards": 0,
    }
    for task in manifest["tasks"]:
        checkpoint_path = Path(task["output_root"]) / "checkpoint.json"
        if not checkpoint_path.is_file():
            totals["pending_shards"] += 1
            continue
        checkpoint = read_json(checkpoint_path)
        totals["successful"] += int(checkpoint.get("successful", 0))
        totals["failed_molecules"] += int(checkpoint.get("failed", 0))
        state = checkpoint.get("status")
        if state == "DONE":
            totals["completed_shards"] += 1
        elif state == "RUNNING":
            totals["running_shards"] += 1
        elif state in {"FAILED", "INTERRUPTED"}:
            task_index = int(task["task_index"])
            attempts = [
                read_json(_attempt_path(manifest_path, task_index, attempt))
                for attempt in range(1, MAX_ATTEMPTS + 1)
                if _attempt_path(manifest_path, task_index, attempt).is_file()
            ]
            if attempts and (attempts[-1].get("retryable") is False or attempts[-1]["attempt"] >= MAX_ATTEMPTS):
                totals["failed_shards"] += 1
            else:
                totals["pending_shards"] += 1
        else:
            totals["pending_shards"] += 1
    processed = totals["successful"] + totals["failed_molecules"]
    remaining = max(0, totals["expected"] - processed)
    percentage = 100.0 * processed / totals["expected"] if totals["expected"] else 100.0
    text = "\n".join([
        f"timestamp={utc_now()}",
        f"campaign={manifest['campaign']}",
        f"total_expected={totals['expected']}",
        f"processed={processed}",
        f"successful={totals['successful']}",
        f"failed_molecules={totals['failed_molecules']}",
        f"remaining={remaining}",
        f"percentage={percentage:.3f}",
        f"completed_shards={totals['completed_shards']}",
        f"running_shards={totals['running_shards']}",
        f"pending_shards={totals['pending_shards']}",
        f"failed_shards={totals['failed_shards']}",
    ]) + "\n"
    path = manifest_path.parent / "production_progress.txt"
    atomic_write_text(path, text)
    return text, path


def emit_row(manifest: dict[str, Any], task_index: int) -> None:
    task = manifest["tasks"][task_index]
    scientific = manifest["scientific"]
    resources = manifest["resources"]
    slurm = manifest["slurm"]
    fields = [
        task["input_path"], task["output_root"], manifest["scratch_root"], scientific["mode"],
        scientific["label_column"], resources["batch_size"], resources["prep_workers"],
        resources["java_threads"], resources["java_heap"], resources["java_lifecycle"],
        manifest["max_attempts"], slurm["conda_root"], slurm["conda_env"],
        1 if slurm["staging_enabled"] else 0, manifest["campaign"],
        1 if resources["retain_scientific_artifacts"] else 0, manifest["project_root"],
        manifest.get("io", {}).get("input_format", "csv"),
        manifest.get("io", {}).get("input_table") or "",
        manifest.get("io", {}).get("input_query") or "",
        manifest.get("io", {}).get("id_column", "MOLECULE_NAME"),
        manifest.get("io", {}).get("smiles_column", "SMILES"),
        manifest.get("io", {}).get("output_format", "csv"),
        manifest.get("io", {}).get("output_table", "demiurge_features"),
        manifest.get("io", {}).get("metadata_table", "demiurge_metadata"),
    ]
    sys.stdout.buffer.write(b"\0".join(str(value).encode("utf-8") for value in fields) + b"\0")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    submit = commands.add_parser("submit", help="create a manifest and submit a bounded SLURM retry chain")
    submit.add_argument("--project-root", type=Path, required=True, help="deployed Demiurge source root")
    submit.add_argument("--input-dir", type=Path, required=True, help="directory containing input shards")
    submit.add_argument("--pattern", default="*.csv", help="input filename glob (default: *.csv)")
    submit.add_argument("--input-format", choices=("csv", "sqlite"), default="csv", help="shard format")
    sqlite_source = submit.add_mutually_exclusive_group()
    sqlite_source.add_argument("--input-table", help="SQLite table/view; exclusive with --input-query")
    sqlite_source.add_argument("--input-query", help="read-only SQLite SELECT; exclusive with --input-table")
    submit.add_argument("--id-column", default="MOLECULE_NAME", help="molecule identifier column")
    submit.add_argument("--smiles-column", default="SMILES", help="SMILES column")
    submit.add_argument("--output-format", choices=("csv", "sqlite"), default="csv", help="one output per shard")
    submit.add_argument("--output-table", default="demiurge_features", help="SQLite result table")
    submit.add_argument("--metadata-table", default="demiurge_metadata", help="SQLite metadata table")
    submit.add_argument("--output-root", type=Path, required=True, help="durable campaign parent")
    submit.add_argument("--scratch-root", type=Path, required=True, help="compute-node staging/scratch parent")
    submit.add_argument("--campaign", required=True, help="new campaign directory name")
    submit.add_argument("--mode", choices=("1H", "13C", "FP", "hybrid", "total"), help="feature mode")
    submit.add_argument("--label-column", type=column_selector, help="one-based position or column name")
    submit.add_argument("--batch-size", type=positive_integer, help="molecules per atomic batch")
    submit.add_argument("--prep-workers", type=positive_integer, help="parallel preparation processes")
    submit.add_argument("--java-threads", type=positive_integer, help="predictor threads per nucleus")
    submit.add_argument("--java-heap", help="Java maximum heap, for example 4G")
    submit.add_argument("--java-lifecycle", choices=("persistent", "per-batch"), help="predictor lifecycle")
    submit.add_argument("--cpus-per-task", type=positive_integer, help="SLURM CPUs per worker")
    submit.add_argument("--memory", help="SLURM memory per worker")
    submit.add_argument("--partition", help="SLURM partition")
    submit.add_argument("--time", help="SLURM wall-time limit")
    submit.add_argument("--max-concurrent", type=positive_integer, help="array concurrency cap")
    submit.add_argument("--job-name-prefix", help="sanitized SLURM job-name prefix")
    submit.add_argument("--conda-root", help="Conda installation root")
    submit.add_argument("--conda-env", help="Conda environment name")
    submit.add_argument("--no-staging", action="store_true", help="read directly instead of staging to scratch")
    submit.add_argument("--retain-scientific-artifacts", action="store_true", help="retain MOL and raw NMR files")
    submit.add_argument("--dry-run", action="store_true", help="write manifest and commands without sbatch")
    resume = commands.add_parser("resume", help="submit remaining compatible shards")
    resume.add_argument("--manifest", type=Path, required=True)
    resume.add_argument("--dry-run", action="store_true")
    status = commands.add_parser("status", help="read durable campaign progress")
    status.add_argument("--manifest", type=Path, required=True)
    row = commands.add_parser("row")
    row.add_argument("--manifest", type=Path, required=True)
    row.add_argument("--task-index", type=int, required=True)
    row.add_argument("--attempt", type=int, required=True)
    decision = commands.add_parser("decision")
    decision.add_argument("--manifest", type=Path, required=True)
    decision.add_argument("--task-index", type=int, required=True)
    decision.add_argument("--attempt", type=int, required=True)
    record = commands.add_parser("record")
    record.add_argument("--manifest", type=Path, required=True)
    record.add_argument("--task-index", type=int, required=True)
    record.add_argument("--attempt", type=int, required=True)
    record.add_argument("--exit-code", type=int, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "submit":
        path, manifest = prepare_manifest(args)
        jobs = submit_chain(path, manifest, list(range(len(manifest["tasks"]))), dry_run=args.dry_run)
        print(f"manifest={path}")
        print("job_ids=" + ",".join(jobs))
        return 0
    path, manifest = load_manifest(args.manifest)
    if args.command == "row":
        emit_row(manifest, args.task_index)
    elif args.command == "decision":
        print(retry_decision(path, manifest, args.task_index, args.attempt))
    elif args.command == "record":
        print(json.dumps(record_attempt(path, manifest, args.task_index, args.attempt, args.exit_code), sort_keys=True))
    elif args.command == "status":
        text, output = campaign_status(path, manifest)
        print(text, end="")
        print(f"progress_file={output}")
    else:
        ensure_previous_arrays_inactive(manifest)
        completed = {
            int(task["task_index"]) for task in manifest["tasks"]
            if (Path(task["output_root"]) / "checkpoint.json").is_file()
            and read_json(Path(task["output_root"]) / "checkpoint.json").get("status") == "DONE"
        }
        remaining = [int(task["task_index"]) for task in manifest["tasks"] if int(task["task_index"]) not in completed]
        used = max((int(item["attempt"]) for item in manifest.get("submissions", [])), default=0)
        if not remaining:
            print("No incomplete shards")
        elif used >= MAX_ATTEMPTS:
            raise RuntimeError("Maximum campaign attempts have already been submitted")
        else:
            jobs = submit_chain(path, manifest, remaining, start_attempt=used + 1, dry_run=args.dry_run)
            print("job_ids=" + ",".join(jobs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
