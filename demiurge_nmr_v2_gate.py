#!/usr/bin/env python3
"""Zero-tolerance cross-repository SPECTRAPRINTS_NMR_V2 parity gate."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

from demiurge_bin.bucketing import nmr_vector, parse_prediction_csv, validate_prediction_records
from demiurge_bin.preparation import prepare_mol_v2 as demiurge_prepare
from demiurge_bin.run_state import atomic_write_json, file_sha256


BRANCH_FILE = ".spectraprints_unified_profile_branches.jsonl"


def _load_corpus(path: Path) -> list[dict[str, str]]:
    rows = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        value = json.loads(line)
        molecule_id = value.get("molecule_id") or value.get("internal_id")
        if not molecule_id or "smiles" not in value:
            raise ValueError(f"Invalid parity corpus row {line_number}")
        rows.append({"molecule_id": str(molecule_id), "smiles": str(value["smiles"])})
    if not rows:
        raise ValueError("Parity corpus is empty")
    if len({row["molecule_id"] for row in rows}) != len(rows):
        raise ValueError("Parity corpus contains duplicate molecule IDs")
    return rows


def _screen_modules(screen_root: Path):
    root = screen_root.expanduser().resolve()
    if not (root / "engine" / "training_contract.py").is_file():
        raise FileNotFoundError(f"Not a screen_SPECTRAprints checkout: {root}")
    sys.path.insert(0, str(root))
    try:
        training = importlib.import_module("engine.training_contract")
        java = importlib.import_module("engine.demiurge_bin.predictor")
        parser = importlib.import_module("engine.demiurge_sqlite_batch")
    finally:
        sys.path.pop(0)
    return training, java, parser


def _write_mols(
    rows: list[dict[str, str]],
    output: Path,
    prepare: Callable[[str], tuple[str, str]],
) -> tuple[dict[str, str], list[dict[str, str]]]:
    mols = output / "mols"
    mols.mkdir()
    canonical: dict[str, str] = {}
    failures: list[dict[str, str]] = []
    for row in rows:
        try:
            identity, block = prepare(row["smiles"])
            canonical[row["molecule_id"]] = identity
            with (mols / f"{row['molecule_id']}.mol").open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(block)
        except Exception as exc:
            failures.append({
                "molecule_id": row["molecule_id"],
                "error_type": type(exc).__name__,
                "error_message": str(exc),
            })
    atomic_write_json(output / "canonical_identity.json", canonical)
    atomic_write_json(output / "preparation_failures.json", {"failures": failures})
    return canonical, failures


def _emit_vectors(
    molecule_ids: list[str],
    raw_h: Path,
    raw_c: Path,
    *,
    implementation: str,
    screen_training: Any = None,
    screen_parser: Any = None,
) -> None:
    vectors: dict[str, Any] = {}
    failures: list[dict[str, str]] = []
    for molecule_id in molecule_ids:
        try:
            h_path = raw_h / f"{molecule_id}.csv"
            c_path = raw_c / f"{molecule_id}.csv"
            if implementation == "screen":
                h_records = screen_parser.parse_prediction_csv(h_path, "1H")
                c_records = screen_parser.parse_prediction_csv(c_path, "13C")
                h_shifts = validate_prediction_records(h_records, "1H")
                c_shifts = validate_prediction_records(c_records, "13C")
                combined, diagnostics = screen_training.nmr_vector(h_shifts, c_shifts)
            else:
                h_shifts = validate_prediction_records(parse_prediction_csv(h_path, "1H"), "1H")
                c_shifts = validate_prediction_records(parse_prediction_csv(c_path, "13C"), "13C")
                combined, diagnostics = nmr_vector(h_shifts, c_shifts)
            diagnostics = {
                "h_in_range": int(diagnostics["h_in_range"]),
                "h_out_of_range": int(diagnostics["h_out_of_range"]),
                "c_in_range": int(diagnostics["c_in_range"]),
                "c_out_of_range": int(diagnostics["c_out_of_range"]),
            }
            vector = [int(value) for value in combined.tolist()]
            vectors[molecule_id] = {
                "1H": vector[:200],
                "13C": vector[200:],
                "H_C": vector,
                "diagnostics": diagnostics,
            }
        except Exception as exc:
            failures.append({
                "molecule_id": molecule_id,
                "error_type": type(exc).__name__,
                "error_message": str(exc),
            })
    atomic_write_json(raw_h.parent / "nmr_vectors.json", vectors)
    atomic_write_json(raw_h.parent / "prediction_failures.json", {"failures": failures})


def emit(args: argparse.Namespace) -> int:
    if args.implementation == "screen" and args.screen_root is None:
        raise ValueError("--screen-root is required for --implementation screen")
    output = args.output_root.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"Parity output root must be empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    rows = _load_corpus(args.corpus)
    screen_training = screen_java = screen_parser = None
    if args.implementation == "screen":
        screen_training, screen_java, screen_parser = _screen_modules(args.screen_root)
        prepare = screen_training.prepare_mol_v2
        java = screen_java
    else:
        from demiurge_bin import predictor as java
        prepare = demiurge_prepare

    canonical, prep_failures = _write_mols(rows, output, prepare)
    successful_ids = [row["molecule_id"] for row in rows if row["molecule_id"] in canonical]
    os.environ["SPECTRAPRINTS_JAVA_LIFECYCLE"] = args.java_lifecycle
    os.environ["SPECTRAPRINTS_JAVA_PREDICTOR_MODE"] = "thread-local"
    os.environ["SPECTRAPRINTS_UNIFIED_PROFILE"] = "1"
    os.environ["SPECTRAPRINTS_JAVA_DIAGNOSTICS_DIR"] = str(output / "diagnostics")
    raw_h = output / "raw_1h"
    raw_c = output / "raw_13c"
    raw_h.mkdir()
    raw_c.mkdir()
    try:
        if successful_ids:
            if args.implementation == "screen":
                old_cwd = Path.cwd()
                os.chdir(output)
                try:
                    h_result = java.run_java_batch_processor(output / "mols", "1H", args.java_threads, args.java_heap)
                    c_result = java.run_java_batch_processor(output / "mols", "13C", args.java_threads, args.java_heap)
                finally:
                    os.chdir(old_cwd)
                if h_result is None or c_result is None:
                    raise RuntimeError("Screen Java predictor failed")
                generated_h = output / "predicted_spectra_1H"
                generated_c = output / "predicted_spectra_13C"
                raw_h.rmdir(); raw_c.rmdir()
                generated_h.rename(raw_h); generated_c.rename(raw_c)
            else:
                h_result = java.run_java_batch_processor(output / "mols", "1H", args.java_threads, args.java_heap, output_directory=raw_h)
                c_result = java.run_java_batch_processor(output / "mols", "13C", args.java_threads, args.java_heap, output_directory=raw_c)
                if h_result is None or c_result is None:
                    raise RuntimeError("Demiurge Java predictor failed")
    finally:
        java.shutdown_persistent_java_processors("parity-gate-finally")
    _emit_vectors(
        successful_ids,
        raw_h,
        raw_c,
        implementation=args.implementation,
        screen_training=screen_training,
        screen_parser=screen_parser,
    )
    atomic_write_json(output / "gate_manifest.json", {
        "schema_version": 1,
        "implementation": args.implementation,
        "corpus_sha256": file_sha256(args.corpus),
        "molecules": len(rows),
        "prepared": len(successful_ids),
        "preparation_failed": len(prep_failures),
        "java_threads": args.java_threads,
        "java_heap": args.java_heap,
        "java_lifecycle": args.java_lifecycle,
    })
    return 0


def _compare_file(left: Path, right: Path, label: str) -> None:
    if not left.is_file() or not right.is_file():
        raise RuntimeError(f"Missing {label}: {left} or {right}")
    if left.read_bytes() != right.read_bytes():
        raise RuntimeError(f"Exact {label} mismatch: {left} vs {right}")


def _compare_tree(left: Path, right: Path, label: str, pattern: str = "*") -> None:
    left_files = sorted(path.relative_to(left) for path in left.rglob(pattern) if path.is_file())
    right_files = sorted(path.relative_to(right) for path in right.rglob(pattern) if path.is_file())
    if left_files != right_files:
        raise RuntimeError(f"{label} filename set differs: {left_files} vs {right_files}")
    for relative in left_files:
        _compare_file(left / relative, right / relative, f"{label}/{relative}")


def compare(args: argparse.Namespace) -> int:
    left = args.screen_output.expanduser().resolve()
    right = args.demiurge_output.expanduser().resolve()
    for name in ("canonical_identity.json", "preparation_failures.json", "nmr_vectors.json", "prediction_failures.json"):
        left_value = json.loads((left / name).read_text(encoding="utf-8"))
        right_value = json.loads((right / name).read_text(encoding="utf-8"))
        if left_value != right_value:
            raise RuntimeError(f"Exact scientific JSON mismatch: {name}")
    _compare_tree(left / "mols", right / "mols", "V3000 MOL", "*.mol")
    _compare_tree(left / "raw_1h", right / "raw_1h", "raw 1H", "*.csv")
    _compare_tree(left / "raw_13c", right / "raw_13c", "raw 13C", "*.csv")
    for nucleus in ("raw_1h", "raw_13c"):
        left_branch = left / nucleus / BRANCH_FILE
        right_branch = right / nucleus / BRANCH_FILE
        _compare_file(left_branch, right_branch, f"{nucleus} 3D branch status")
    print("PASS exact SPECTRAPRINTS_NMR_V2 parity")
    return 0


def _collect_run_artifacts(root: Path, category: str) -> dict[str, bytes]:
    collected: dict[str, bytes] = {}
    for batch in sorted((root / "batches").glob("batch_*")):
        directory = batch / "scientific_artifacts" / category
        if not directory.is_dir():
            raise RuntimeError(f"Run did not retain required {category} artifacts: {directory}")
        for path in directory.glob("*"):
            if path.is_file() and not path.name.startswith("."):
                if path.name in collected:
                    raise RuntimeError(f"Duplicate retained artifact {category}/{path.name}")
                collected[path.name] = path.read_bytes()
    return collected


def _collect_run_branches(root: Path, category: str) -> dict[str, list[dict[str, Any]]]:
    collected: dict[str, list[dict[str, Any]]] = {}
    for batch in sorted((root / "batches").glob("batch_*")):
        path = batch / "scientific_artifacts" / category / BRANCH_FILE
        if not path.is_file():
            raise RuntimeError(f"Run did not retain required 3D branch status: {path}")
        values = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
        key = batch.name
        collected[key] = sorted(values, key=lambda value: (value["internal_id"], value["nucleus"]))
    return collected


def _collect_run_metadata(root: Path) -> list[dict[str, Any]]:
    values = []
    for batch in sorted((root / "batches").glob("batch_*")):
        for line in (batch / "metadata.jsonl").read_text(encoding="utf-8").splitlines():
            values.append(json.loads(line))
    return sorted(values, key=lambda item: int(item["source_index"]))


def compare_runs(args: argparse.Namespace) -> int:
    left = args.left_output.expanduser().resolve()
    right = args.right_output.expanduser().resolve()
    left_summary = json.loads((left / "summary.json").read_text(encoding="utf-8"))
    right_summary = json.loads((right / "summary.json").read_text(encoding="utf-8"))
    if left_summary["contract_id"] != right_summary["contract_id"]:
        raise RuntimeError("Run scientific contract IDs differ")
    for key in ("total", "successful", "failed"):
        if left_summary[key] != right_summary[key]:
            raise RuntimeError(f"Run summary scientific count differs: {key}")
    _compare_file(Path(left_summary["final_output"]), Path(right_summary["final_output"]), "final feature CSV")
    _compare_file(left / "failures.jsonl", right / "failures.jsonl", "failure identities/reasons")
    if _collect_run_metadata(left) != _collect_run_metadata(right):
        raise RuntimeError("Run canonical identities, feature hashes or failure metadata differ")
    for category in ("mols", "raw_1h", "raw_13c"):
        left_artifacts = _collect_run_artifacts(left, category)
        right_artifacts = _collect_run_artifacts(right, category)
        if left_artifacts != right_artifacts:
            raise RuntimeError(f"Exact retained run artifacts differ: {category}")
    for category in ("raw_1h", "raw_13c"):
        if _collect_run_branches(left, category) != _collect_run_branches(right, category):
            raise RuntimeError(f"Exact retained 3D branch status differs: {category}")
    print("PASS exact local/SLURM-worker and lifecycle scientific parity")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    emit_parser = commands.add_parser("emit")
    emit_parser.add_argument("--implementation", choices=("screen", "demiurge"), required=True)
    emit_parser.add_argument("--screen-root", type=Path)
    emit_parser.add_argument("--corpus", type=Path, required=True)
    emit_parser.add_argument("--output-root", type=Path, required=True)
    emit_parser.add_argument("--java-threads", type=int, default=2)
    emit_parser.add_argument("--java-heap", default="4G")
    emit_parser.add_argument("--java-lifecycle", choices=("per-batch", "persistent"), default="persistent")
    compare_parser = commands.add_parser("compare")
    compare_parser.add_argument("--screen-output", type=Path, required=True)
    compare_parser.add_argument("--demiurge-output", type=Path, required=True)
    run_parser = commands.add_parser("compare-runs")
    run_parser.add_argument("--left-output", type=Path, required=True)
    run_parser.add_argument("--right-output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "emit":
        if args.java_threads <= 0:
            raise ValueError("java_threads must be positive")
        return emit(args)
    if args.command == "compare":
        return compare(args)
    return compare_runs(args)


if __name__ == "__main__":
    raise SystemExit(main())
