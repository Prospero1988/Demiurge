"""Exact indexed-CSV parsing and NMR V2 count bucketing."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Optional

import numpy as np

from .contracts import C_BINS, C_MAX, C_MIN, H_BINS, H_MAX, H_MIN


PredictionRecord = tuple[int, str, Optional[int], float]


def parse_prediction_csv(csv_path: Path, nucleus: str) -> list[PredictionRecord]:
    if nucleus not in {"1H", "13C"}:
        raise ValueError(f"Unsupported nucleus: {nucleus}")
    required = (
        {"mol_atom_index", "element", "parent_atom_index", "shift"}
        if nucleus == "1H"
        else {"mol_atom_index", "element", "shift"}
    )
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{csv_path} is missing columns: {sorted(missing)}")
        records: list[PredictionRecord] = []
        for row in reader:
            parent_raw = row.get("parent_atom_index")
            parent = None if parent_raw is None or str(parent_raw).strip() == "" else int(parent_raw)
            records.append((int(row["mol_atom_index"]), str(row["element"]), parent, float(row["shift"])))
    return records


def validate_prediction_records(records: list[PredictionRecord], nucleus: str) -> list[float]:
    expected_element = "H" if nucleus == "1H" else "C"
    seen: set[int] = set()
    shifts: list[float] = []
    for atom_index, element, parent_index, shift in records:
        atom = int(atom_index)
        if atom <= 0 or atom in seen:
            raise ValueError(f"Duplicate/non-positive {nucleus} mol_atom_index {atom}")
        seen.add(atom)
        if element != expected_element:
            raise ValueError(f"Unexpected {nucleus} element {element!r}")
        if nucleus == "1H" and parent_index is not None and int(parent_index) <= 0:
            raise ValueError(f"Invalid 1H parent_atom_index {parent_index}")
        value = float(shift)
        if not math.isfinite(value):
            raise ValueError(f"Non-finite {nucleus} shift")
        shifts.append(value)
    return shifts


def bucket_shifts(shifts: list[float], ppm_min: float, ppm_max: float, n_bins: int) -> tuple[np.ndarray, int, int]:
    if n_bins <= 0 or ppm_max <= ppm_min:
        raise ValueError("Invalid NMR bucket contract")
    buckets = np.zeros(int(n_bins), dtype=np.float32)
    inside = outside = 0
    width = (float(ppm_max) - float(ppm_min)) / float(n_bins)
    for shift in shifts:
        value = float(shift)
        if value < ppm_min or value > ppm_max:
            outside += 1
            continue
        index = n_bins - 1 if value == ppm_max else int((value - ppm_min) // width)
        buckets[index] += 1.0
        inside += 1
    return buckets, inside, outside


def nmr_vector(h_shifts: list[float], c_shifts: list[float]) -> tuple[np.ndarray, dict[str, int]]:
    h, h_inside, h_outside = bucket_shifts(h_shifts, H_MIN, H_MAX, H_BINS)
    c, c_inside, c_outside = bucket_shifts(c_shifts, C_MIN, C_MAX, C_BINS)
    return np.concatenate([h, c]).astype(np.float32, copy=False), {
        "h_total": len(h_shifts),
        "c_total": len(c_shifts),
        "h_in_range": h_inside,
        "c_in_range": c_inside,
        "h_out_of_range": h_outside,
        "c_out_of_range": c_outside,
    }


def vectors_from_raw(h_path: Path, c_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    h_shifts = validate_prediction_records(parse_prediction_csv(h_path, "1H"), "1H")
    c_shifts = validate_prediction_records(parse_prediction_csv(c_path, "13C"), "13C")
    combined, diagnostics = nmr_vector(h_shifts, c_shifts)
    return combined[:H_BINS], combined[H_BINS:], combined, diagnostics
