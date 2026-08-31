"""Marker-owned input and workspace staging for SLURM workers."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import uuid
from pathlib import Path

from demiurge_bin.run_state import atomic_write_json, file_sha256, utc_now


OWNER_FILE = ".demiurge-staging-owner.json"


def _safe_child(root: Path, child: Path) -> Path:
    resolved_root = root.expanduser().resolve()
    if resolved_root == Path(resolved_root.anchor):
        raise RuntimeError(f"Refusing filesystem root as staging root: {resolved_root}")
    resolved = child.resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise RuntimeError(f"Staging directory escapes root: {resolved}") from exc
    return resolved


def stage_input(
    source: Path,
    scratch_root: Path,
    campaign: str,
    job_id: str,
    task_index: int,
    attempt: int,
) -> dict[str, str]:
    source = source.expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Staging source does not exist: {source}")
    root = scratch_root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    safe_campaign = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in campaign)
    token = uuid.uuid4().hex
    directory = root / f"demiurge_{safe_campaign}_j{job_id}_t{task_index:04d}_a{attempt:02d}_{token[:8]}"
    _safe_child(root, directory)
    directory.mkdir(mode=0o700)
    marker = {
        "schema_version": 1,
        "owner_token": token,
        "source": str(source),
        "source_sha256": file_sha256(source),
        "created_at": utc_now(),
    }
    atomic_write_json(directory / OWNER_FILE, marker)
    runtime_input = directory / source.name
    try:
        shutil.copy2(source, runtime_input)
        if file_sha256(runtime_input) != marker["source_sha256"]:
            raise RuntimeError("Staged input hash mismatch")
    except Exception:
        cleanup_staging(root, directory, token)
        raise
    return {
        "staging_directory": str(directory),
        "runtime_input": str(runtime_input),
        "owner_token": token,
    }


def cleanup_staging(scratch_root: Path, directory: Path, token: str) -> None:
    root = scratch_root.expanduser().resolve()
    lexical = directory.absolute()
    if lexical.is_symlink():
        raise RuntimeError(f"Refusing symlink staging cleanup: {lexical}")
    resolved = _safe_child(root, lexical)
    if not resolved.exists():
        return
    marker_path = resolved / OWNER_FILE
    if not marker_path.is_file():
        raise RuntimeError(f"Staging owner marker is missing: {marker_path}")
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if marker.get("owner_token") != token:
        raise RuntimeError(f"Staging owner marker mismatch: {marker_path}")
    shutil.rmtree(resolved)


def main() -> int:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    stage = commands.add_parser("stage")
    stage.add_argument("--source", type=Path, required=True)
    stage.add_argument("--scratch-root", type=Path, required=True)
    stage.add_argument("--campaign", required=True)
    stage.add_argument("--job-id", required=True)
    stage.add_argument("--task-index", type=int, required=True)
    stage.add_argument("--attempt", type=int, required=True)
    cleanup = commands.add_parser("cleanup")
    cleanup.add_argument("--scratch-root", type=Path, required=True)
    cleanup.add_argument("--directory", type=Path, required=True)
    cleanup.add_argument("--owner-token", required=True)
    args = parser.parse_args()
    if args.command == "stage":
        result = stage_input(args.source, args.scratch_root, args.campaign, args.job_id, args.task_index, args.attempt)
        print(result["staging_directory"])
        print(result["runtime_input"])
        print(result["owner_token"])
    else:
        cleanup_staging(args.scratch_root, args.directory, args.owner_token)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
