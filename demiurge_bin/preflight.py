"""Fail-closed, inexpensive deployment preflight for Demiurge."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from .contracts import verify_predictor_artifacts
from .predictor import _get_build_dir, _get_project_root, _resolve_java_tool


def _tool_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "-version"],
        check=True,
        capture_output=True,
        text=True,
    )
    output = (completed.stdout + "\n" + completed.stderr).strip()
    if not output:
        raise RuntimeError(f"Java tool produced no version information: {executable}")
    return output.splitlines()[0]


def run_preflight(project_root: Path | None = None) -> dict[str, Any]:
    """Validate the complete standalone Java toolchain before expensive work."""
    root = (project_root or _get_project_root()).expanduser().resolve()
    java = _resolve_java_tool("java")
    javac = _resolve_java_tool("javac")
    artifacts = verify_predictor_artifacts(root)
    build = _get_build_dir()
    return {
        "status": "PASS",
        "project_root": str(root),
        "java": {"path": java, "version": _tool_version(java)},
        "javac": {"path": javac, "version": _tool_version(javac)},
        "predictor_artifact_sha256": artifacts,
        "java_build_directory": str(build),
        "java_build_directory_writable": True,
    }
