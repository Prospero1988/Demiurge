"""Backward-compatible streaming CSV input and output."""

from __future__ import annotations

import csv
import json
import os
import shutil
from pathlib import Path
from typing import Any, Iterator

import pandas as pd

from demiurge_bin.run_state import atomic_write_json, atomic_write_text, file_sha256

from .base import ColumnSelector, InputDescription, InputReader, OutputWriter, record_from_values, resolve_column


class CsvInputReader(InputReader):
    def __init__(self, path: Path, id_column: str, smiles_column: str, label_column: ColumnSelector):
        self.path = path
        self.id_column = id_column
        self.smiles_column = smiles_column
        self.label_column = label_column
        self._description: InputDescription | None = None
        self._label_index: int | None = None

    @staticmethod
    def _read(path: Path, **kwargs: Any):
        return pd.read_csv(path, sep=None, engine="python", **kwargs)

    def describe(self) -> InputDescription:
        if self._description is not None:
            return self._description
        header = self._read(self.path, nrows=0)
        columns = [str(column) for column in header.columns]
        for option, name in (("id_column", self.id_column), ("smiles_column", self.smiles_column)):
            if columns.count(name) != 1:
                reason = "does not exist" if name not in columns else "is ambiguous"
                raise ValueError(f"{option}={name!r} {reason} in the input")
        self._label_index, label_name = resolve_column(columns, self.label_column, "label_column")
        total = 0
        with self._read(self.path, usecols=[self.id_column], chunksize=100_000) as chunks:
            for chunk in chunks:
                total += len(chunk)
        self._description = InputDescription(
            total_rows=total,
            label_name=label_name,
            configuration={
                "input_format": "csv",
                "id_column": self.id_column,
                "smiles_column": self.smiles_column,
                "label_column": self.label_column,
                "label_name": label_name,
            },
        )
        return self._description

    def iter_batches(self, batch_size: int, start_index: int = 0) -> Iterator[tuple[int, list[dict[str, Any]]]]:
        self.describe()
        assert self._label_index is not None
        pending: list[dict[str, Any]] = []
        pending_offset: int | None = None
        source_index = 0
        with self._read(self.path, chunksize=batch_size) as chunks:
            for frame in chunks:
                for _, row in frame.iterrows():
                    current = source_index
                    source_index += 1
                    if current < start_index:
                        continue
                    if pending_offset is None:
                        pending_offset = current
                    pending.append(record_from_values(
                        current,
                        row.get(self.id_column),
                        row.get(self.smiles_column),
                        row.iloc[self._label_index],
                        self.describe().label_name,
                        self.id_column,
                    ))
                    if len(pending) == batch_size:
                        yield pending_offset, pending
                        pending = []
                        pending_offset = None
        if pending:
            assert pending_offset is not None
            yield pending_offset, pending


class CsvOutputWriter(OutputWriter):
    def __init__(self, output_root: Path, input_stem: str, mode: str, feature_count: int):
        self.output_root = output_root
        self.input_stem = input_stem
        self.mode = mode
        self.feature_count = feature_count

    def initialize(self, *, resume: bool) -> None:
        del resume

    def configuration(self) -> dict[str, Any]:
        return {"output_format": "csv"}

    def commit_batch(
        self,
        batch_index: int,
        feature_rows: list[dict[str, Any]],
        metadata: list[dict[str, Any]],
        timings: dict[str, float],
        scratch_batch: Path,
        retain_scientific_artifacts: bool,
    ) -> Path:
        batches = self.output_root / "batches"
        batches.mkdir(parents=True, exist_ok=True)
        final = batches / f"batch_{batch_index:08d}"
        temporary = batches / f".batch_{batch_index:08d}.{os.getpid()}.tmp"
        if temporary.exists():
            shutil.rmtree(temporary)
        temporary.mkdir()
        feature_path = temporary / "features.csv"
        with feature_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["MOLECULE_NAME", "LABEL"] + [
                f"FEATURE_{index}" for index in range(1, self.feature_count + 1)
            ])
            for row in feature_rows:
                writer.writerow([row["molecule_name"], row["label"], *row["vector"]])
        metadata_path = temporary / "metadata.jsonl"
        atomic_write_text(metadata_path, "".join(
            json.dumps(item, sort_keys=True, ensure_ascii=False) + "\n" for item in metadata
        ))
        _retain_artifacts(temporary, scratch_batch, retain_scientific_artifacts)
        commit = {
            "schema_version": 1,
            "batch_index": batch_index,
            "successful": len(feature_rows),
            "failed": sum(1 for item in metadata if item["status"] == "FAILED"),
            "features_sha256": file_sha256(feature_path),
            "metadata_sha256": file_sha256(metadata_path),
            "timing_seconds": timings,
            "committed_at": _utc_now(),
        }
        atomic_write_json(temporary / "commit.json", commit)
        if final.exists():
            _validate_replay(final, temporary, retain_scientific_artifacts, include_features=True)
            shutil.rmtree(temporary)
            return final
        os.replace(temporary, final)
        return final

    def finalize(self, *, total: int, successful: int, failed: int) -> tuple[Path, Path]:
        del total, successful, failed
        final_directory = self.output_root / "generated_ML_inputs"
        final_directory.mkdir(parents=True, exist_ok=True)
        final = final_directory / f"{self.input_stem}_{self.mode}_ML_input.csv"
        temporary = final.with_name(f".{final.name}.{os.getpid()}.tmp")
        wrote_header = False
        with temporary.open("w", encoding="utf-8", newline="") as destination:
            for batch in sorted((self.output_root / "batches").glob("batch_*")):
                source = batch / "features.csv"
                with source.open(encoding="utf-8", newline="") as handle:
                    for line_index, line in enumerate(handle):
                        if line_index == 0 and wrote_header:
                            continue
                        destination.write(line)
                        wrote_header = True
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(temporary, final)
        return final, _aggregate_failures(self.output_root)


def _utc_now() -> str:
    from demiurge_bin.run_state import utc_now

    return utc_now()


def _retain_artifacts(destination: Path, scratch_batch: Path, enabled: bool) -> None:
    if not enabled:
        return
    artifacts = destination / "scientific_artifacts"
    artifacts.mkdir()
    for name in ("mols", "raw_1h", "raw_13c"):
        source = scratch_batch / name
        if source.is_dir():
            shutil.copytree(source, artifacts / name)


def _validate_replay(final: Path, temporary: Path, retain_artifacts: bool, *, include_features: bool) -> None:
    relatives = [Path("metadata.jsonl")]
    if include_features:
        relatives.insert(0, Path("features.csv"))
    for relative in relatives:
        if (final / relative).read_bytes() != (temporary / relative).read_bytes():
            raise RuntimeError(f"Existing atomic batch differs during replay: {final / relative}")
    for artifact_root in (final / "scientific_artifacts", temporary / "scientific_artifacts"):
        if artifact_root.exists() != retain_artifacts:
            raise RuntimeError(f"Existing batch artifact retention differs during replay: {final}")
    if retain_artifacts:
        left = {path.relative_to(final / "scientific_artifacts"): path.read_bytes() for path in (final / "scientific_artifacts").rglob("*") if path.is_file()}
        right = {path.relative_to(temporary / "scientific_artifacts"): path.read_bytes() for path in (temporary / "scientific_artifacts").rglob("*") if path.is_file()}
        if left != right:
            raise RuntimeError(f"Existing atomic scientific artifacts differ during replay: {final}")


def _aggregate_failures(output_root: Path) -> Path:
    failures = output_root / "failures.jsonl"
    lines: list[str] = []
    for batch in sorted((output_root / "batches").glob("batch_*")):
        for line in (batch / "metadata.jsonl").read_text(encoding="utf-8").splitlines():
            value = json.loads(line)
            if value.get("status") == "FAILED":
                lines.append(json.dumps(value, sort_keys=True, ensure_ascii=False) + "\n")
    atomic_write_text(failures, "".join(lines))
    return failures
