"""Durable, atomic run state shared by local and SLURM worker execution."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CHECKPOINT_SCHEMA_VERSION = 1


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def input_identity(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": stat.st_size,
        "sha256": file_sha256(resolved),
    }


def compatible_identity(expected: dict[str, Any], current: dict[str, Any]) -> bool:
    return (
        int(expected.get("size", -1)) == int(current.get("size", -2))
        and str(expected.get("sha256")) == str(current.get("sha256"))
    )


def initial_checkpoint(
    *,
    run_id: str,
    input_info: dict[str, Any],
    scientific_config: dict[str, Any],
    total_rows: int,
    io_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    now = utc_now()
    checkpoint = {
        "checkpoint_schema_version": CHECKPOINT_SCHEMA_VERSION,
        "run_id": run_id,
        "status": "INITIALIZED",
        "input_identity": input_info,
        "scientific_config": scientific_config,
        "total_rows": int(total_rows),
        "next_row_index": 0,
        "committed_batches": 0,
        "successful": 0,
        "failed": 0,
        "attempt": 1,
        "started_at": now,
        "heartbeat": now,
        "finished_at": None,
        "failure_type": None,
        "failure_message": None,
    }
    if io_config is not None:
        checkpoint["io_config"] = io_config
    return checkpoint


def validate_resume(
    checkpoint: dict[str, Any],
    *,
    input_info: dict[str, Any],
    scientific_config: dict[str, Any],
    io_config: dict[str, Any] | None = None,
) -> None:
    if int(checkpoint.get("checkpoint_schema_version", -1)) != CHECKPOINT_SCHEMA_VERSION:
        raise RuntimeError("Checkpoint schema is incompatible")
    if not compatible_identity(checkpoint.get("input_identity") or {}, input_info):
        raise RuntimeError("Input content identity differs from checkpoint")
    if checkpoint.get("scientific_config") != scientific_config:
        raise RuntimeError("Scientific configuration differs from checkpoint")
    if io_config is not None and checkpoint.get("io_config") not in (None, io_config):
        raise RuntimeError("Input/output configuration differs from checkpoint")
    if checkpoint.get("status") == "DONE":
        raise RuntimeError("Run is already DONE")


def progress_text(checkpoint: dict[str, Any]) -> str:
    total = int(checkpoint.get("total_rows", 0))
    successful = int(checkpoint.get("successful", 0))
    failed = int(checkpoint.get("failed", 0))
    processed = successful + failed
    remaining = max(0, total - processed)
    percentage = (100.0 * processed / total) if total else 100.0
    lines = [
        f"timestamp={utc_now()}",
        f"run_id={checkpoint.get('run_id', 'unknown')}",
        f"status={checkpoint.get('status', 'UNKNOWN')}",
        f"total={total}",
        f"processed={processed}",
        f"successful={successful}",
        f"failed={failed}",
        f"remaining={remaining}",
        f"percentage={percentage:.3f}",
        f"committed_batches={int(checkpoint.get('committed_batches', 0))}",
        f"heartbeat={checkpoint.get('heartbeat')}",
    ]
    return "\n".join(lines) + "\n"


def write_progress(output_root: Path, checkpoint: dict[str, Any]) -> Path:
    path = output_root / "production_progress.txt"
    atomic_write_text(path, progress_text(checkpoint))
    return path
