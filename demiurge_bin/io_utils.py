"""Owned scratch lifecycle helpers."""

from __future__ import annotations

import json
import shutil
import uuid
from pathlib import Path

from .run_state import atomic_write_json, utc_now


OWNER_FILE = ".demiurge-scratch-owner.json"


def create_owned_scratch(temp_root: Path, run_id: str) -> tuple[Path, str]:
    root = temp_root.expanduser().resolve()
    if root == Path(root.anchor):
        raise RuntimeError(f"Refusing filesystem root as temp root: {root}")
    root.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    directory = root / f"demiurge_{run_id}_{token[:12]}"
    directory.mkdir(mode=0o700)
    atomic_write_json(directory / OWNER_FILE, {
        "schema_version": 1,
        "owner_token": token,
        "run_id": run_id,
        "created_at": utc_now(),
    })
    return directory, token


def cleanup_owned_scratch(temp_root: Path, directory: Path, token: str) -> None:
    root = temp_root.expanduser().resolve()
    lexical = directory.absolute()
    if lexical.is_symlink():
        raise RuntimeError(f"Refusing to clean symlink scratch directory: {lexical}")
    resolved = lexical.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise RuntimeError(f"Scratch directory escapes configured root: {resolved}") from exc
    marker = resolved / OWNER_FILE
    if not marker.is_file():
        raise RuntimeError(f"Owned scratch marker is missing: {marker}")
    document = json.loads(marker.read_text(encoding="utf-8"))
    if document.get("owner_token") != token:
        raise RuntimeError(f"Owned scratch marker does not match: {marker}")
    shutil.rmtree(resolved)
