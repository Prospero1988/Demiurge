"""Interfaces and shared record handling for Demiurge I/O backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import pandas as pd


ColumnSelector = int | str


@dataclass(frozen=True)
class InputDescription:
    total_rows: int
    label_name: str
    configuration: dict[str, Any]


def resolve_column(columns: list[str], selector: ColumnSelector, option: str) -> tuple[int, str]:
    if isinstance(selector, bool):
        raise ValueError(f"{option} must be a one-based position or column name")
    if isinstance(selector, int):
        index = selector - 1
        if index < 0 or index >= len(columns):
            raise ValueError(f"{option}={selector} is outside the input columns")
        return index, str(columns[index])
    name = str(selector)
    matches = [index for index, column in enumerate(columns) if column == name]
    if len(matches) != 1:
        reason = "does not exist" if not matches else "is ambiguous"
        raise ValueError(f"{option}={name!r} {reason} in the input")
    return matches[0], name


def record_from_values(
    source_index: int,
    molecule_id: Any,
    smiles: Any,
    label: Any,
    label_name: str,
    id_name: str,
) -> dict[str, Any]:
    input_error = None
    if pd.isna(molecule_id) or str(molecule_id).strip() == "":
        input_error = "Missing MOLECULE_NAME" if id_name == "MOLECULE_NAME" else f"Missing molecule identifier in {id_name}"
    elif pd.isna(smiles) or str(smiles).strip() == "":
        input_error = "Missing SMILES"
    elif pd.isna(label):
        input_error = f"Missing label in {label_name}"
    scalar_label = None if pd.isna(label) else label.item() if hasattr(label, "item") else label
    return {
        "source_index": int(source_index),
        "internal_id": f"m{int(source_index):08d}",
        "molecule_name": "" if pd.isna(molecule_id) else str(molecule_id),
        "smiles": "" if pd.isna(smiles) else str(smiles),
        "label": scalar_label,
        "input_error": input_error,
    }


class InputReader(ABC):
    @abstractmethod
    def describe(self) -> InputDescription:
        """Validate the source and return stable input metadata."""

    @abstractmethod
    def iter_batches(self, batch_size: int, start_index: int = 0) -> Iterator[tuple[int, list[dict[str, Any]]]]:
        """Yield bounded record batches beginning at source row start_index."""


class OutputWriter(ABC):
    @abstractmethod
    def initialize(self, *, resume: bool) -> None:
        """Create or validate the durable output target."""

    @abstractmethod
    def commit_batch(
        self,
        batch_index: int,
        feature_rows: list[dict[str, Any]],
        metadata: list[dict[str, Any]],
        timings: dict[str, float],
        scratch_batch: Path,
        retain_scientific_artifacts: bool,
    ) -> Path:
        """Commit one idempotent batch and return its durable commit directory."""

    @abstractmethod
    def finalize(self, *, total: int, successful: int, failed: int) -> tuple[Path, Path]:
        """Finalize output and return (primary output, failures JSONL)."""

    @abstractmethod
    def configuration(self) -> dict[str, Any]:
        """Return the resume-relevant output contract."""
