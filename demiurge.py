#!/usr/bin/env python3
"""Demiurge local CLI using the shared NMR V2 scientific core."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from demiurge_bin.java_heap import DEFAULT_JAVA_HEAP, normalize_java_heap
from demiurge_bin.pipeline import RunConfig, read_status, resume_pipeline, run_pipeline


def positive_integer(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def add_operational_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--temp-root", type=Path, default=Path(tempfile.gettempdir()) / "demiurge")
    parser.add_argument("--prep-workers", type=positive_integer, default=4)
    parser.add_argument("--java-threads", type=positive_integer, default=2)
    parser.add_argument("--java-heap", type=normalize_java_heap, default=DEFAULT_JAVA_HEAP)
    parser.add_argument("--batch-size", type=positive_integer, default=500)
    parser.add_argument("--java-lifecycle", choices=("persistent", "per-batch"), default="persistent")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Demiurge SPECTRAPRINTS_NMR_V2 feature generator")
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="run locally or as a scheduler worker")
    run.add_argument("--input", type=Path, required=True)
    run.add_argument("--canonical-input", type=Path, help=argparse.SUPPRESS)
    run.add_argument("--mode", choices=("1H", "13C", "FP", "hybrid", "total"), required=True)
    run.add_argument("--output-root", type=Path, required=True)
    run.add_argument("--label-column", type=positive_integer, default=3)
    run.add_argument("--max-attempts", type=positive_integer, default=3)
    run.add_argument("--retain-scientific-artifacts", action="store_true")
    run.add_argument("--backend", choices=("local", "slurm-worker"), default="local", help=argparse.SUPPRESS)
    add_operational_arguments(run)

    resume = commands.add_parser("resume", help="resume a compatible local run")
    resume.add_argument("--output-root", type=Path, required=True)
    resume.add_argument("--input", type=Path, help=argparse.SUPPRESS)
    resume.add_argument("--temp-root", type=Path)
    resume.add_argument("--prep-workers", type=positive_integer)
    resume.add_argument("--java-threads", type=positive_integer)
    resume.add_argument("--java-heap", type=normalize_java_heap)
    resume.add_argument("--batch-size", type=positive_integer)
    resume.add_argument("--java-lifecycle", choices=("persistent", "per-batch"))
    resume.add_argument("--backend", choices=("local", "slurm-worker"), default="local", help=argparse.SUPPRESS)

    status = commands.add_parser("status", help="read durable progress without modifying scientific state")
    status.add_argument("--output-root", type=Path, required=True)
    return parser


def legacy_arguments(argv: list[str]) -> list[str] | None:
    """Translate the former top-level CLI into the new local run command."""
    if not argv or argv[0] in {"run", "resume", "status", "-h", "--help"}:
        return None
    legacy = argparse.ArgumentParser(add_help=False)
    legacy.add_argument("--csv_path", required=True)
    legacy.add_argument("--predictor", required=True, choices=("1H", "13C", "FP", "hybrid", "total"))
    legacy.add_argument("--label_column", type=positive_integer, required=True)
    legacy.add_argument("--clean", action="store_true")
    values = legacy.parse_args(argv)
    del values.clean
    return [
        "run", "--input", values.csv_path, "--mode", values.predictor,
        "--output-root", str(Path.cwd()), "--label-column", str(values.label_column),
    ]


def main(argv: list[str] | None = None) -> int:
    selected = list(sys.argv[1:] if argv is None else argv)
    translated = legacy_arguments(selected)
    args = build_parser().parse_args(translated if translated is not None else selected)
    if args.command == "status":
        text, path = read_status(args.output_root)
        print(text, end="")
        print(f"progress_file={path}")
        return 0
    if args.command == "resume":
        summary = resume_pipeline(
            args.output_root,
            args.temp_root,
            prep_workers=args.prep_workers,
            java_threads=args.java_threads,
            java_heap=args.java_heap,
            batch_size=args.batch_size,
            java_lifecycle=args.java_lifecycle,
            backend=args.backend,
            input_path=args.input,
        )
    else:
        summary = run_pipeline(RunConfig(
            input_path=args.input,
            mode=args.mode,
            output_root=args.output_root,
            temp_root=args.temp_root,
            label_column=args.label_column,
            prep_workers=args.prep_workers,
            java_threads=args.java_threads,
            java_heap=args.java_heap,
            batch_size=args.batch_size,
            java_lifecycle=args.java_lifecycle,
            max_attempts=args.max_attempts,
            retain_scientific_artifacts=args.retain_scientific_artifacts,
            backend=args.backend,
            canonical_input_path=args.canonical_input,
        ))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
