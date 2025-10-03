#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
concatenator.py
===============

Utility for merging machine-learning input tables into a single dataset.

Supported modes
---------------
- HYBRID (2 datasets): ¹H and ¹³C pseudo-spectra
- TOTAL  (3 datasets): ¹H, ¹³C and ECFP4 fingerprints

Each input must contain two metadata columns:

    • MOLECULE_NAME
    • LABEL                (target value)

All remaining columns are treated as numeric features. Before merging,
feature columns are temporarily prefixed with source-specific prefixes:
    - "H_"   for ¹H
    - "C_"   for ¹³C
    - "FP_"  for fingerprints (ECFP4)

After concatenation, feature columns are renumbered to a generic scheme:

    MOLECULE_NAME, LABEL, FEATURE_1, FEATURE_2, …

Example
-------
    from concatenator import concatenate

    # HYBRID (2 datasets)
    hybrid_df, _ = concatenate(
        ["spectra_1H.csv", "spectra_13C.csv"],
        output_path="hybrid.csv",
    )

    # TOTAL (3 datasets)
    total_df, _ = concatenate(
        ["H.csv", "C.csv", "fp.csv"],
        output_path="total.csv",
    )
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple, Union, List

import pandas as pd

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
META_COLS: List[str] = ["MOLECULE_NAME", "LABEL"]

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def _to_dataframe(data: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    """Return *data* as a fresh ``pandas.DataFrame``."""
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return pd.read_csv(data)


def _ensure_meta(df: pd.DataFrame) -> None:
    """Validate presence of required metadata columns."""
    missing = [c for c in META_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing} in dataset.")


def _split_meta_features(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split *df* into metadata and feature frames."""
    _ensure_meta(df)
    meta_df = df[META_COLS]
    feat_df = df.drop(columns=META_COLS)
    return meta_df, feat_df


def _prefix_features(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Attach *prefix* to every feature column; keep META_COLS unchanged."""
    meta_df, feat_df = _split_meta_features(df)
    feat_df = feat_df.copy()
    feat_df.columns = [f"{prefix}{col}" for col in feat_df.columns]
    return pd.concat([meta_df, feat_df], axis=1)


def _renumber_features(df: pd.DataFrame) -> pd.DataFrame:
    """Rename every feature column to ``FEATURE_n`` while preserving order."""
    feature_cols = [c for c in df.columns if c not in META_COLS]
    mapping = {old: f"FEATURE_{i + 1}" for i, old in enumerate(feature_cols)}
    return df.rename(columns=mapping)


def _merge_on_meta(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    """INNER JOIN two frames on META_COLS without duplicating META_COLS."""
    return pd.merge(left, right, on=META_COLS, how="inner")


# -----------------------------------------------------------------------------
# Core
# -----------------------------------------------------------------------------
def concatenate(
    datasets: Iterable[Union[str, Path, pd.DataFrame]],
    output_path: Union[str, Path, None] = None,
) -> Tuple[pd.DataFrame, Union[str, None]]:
    """
    Concatenate 2 (HYBRID) or 3 (TOTAL) datasets into a unified table.

    Parameters
    ----------
    datasets
        Iterable containing either:
          • exactly 2 elements  -> interpreted as [¹H, ¹³C]
          • exactly 3 elements  -> interpreted as [¹H, ¹³C, FP]
        Each element can be a CSV path or a ``pandas.DataFrame``.
        Order matters (¹H first, then ¹³C, then FP if provided).
    output_path
        Optional path for saving the merged CSV. If *None*, the file is not written.

    Returns
    -------
    (merged_df, saved_path)
        merged_df : ``pandas.DataFrame`` — the combined dataset.
        saved_path : str | None          — path where the CSV was saved, or *None*.
    """
    ds_list = list(datasets)
    n = len(ds_list)

    if n not in (2, 3):
        raise ValueError(
            "Expected 2 datasets (HYBRID: ¹H, ¹³C) or 3 datasets (TOTAL: ¹H, ¹³C, FP). "
            f"Got {n}."
        )

    # Load all as DataFrames
    frames = list(map(_to_dataframe, ds_list))

    # Prefix features according to position:
    # 0 -> ¹H, 1 -> ¹³C, 2 -> FP (if present)
    prefixed = []
    prefixes = ["H_", "C_"] + (["FP_"] if n == 3 else [])
    for df, pref in zip(frames, prefixes):
        prefixed.append(_prefix_features(df, pref))

    # Merge sequentially on META_COLS (INNER)
    merged = prefixed[0]
    for nxt in prefixed[1:]:
        merged = _merge_on_meta(merged, nxt)

    # Renumber feature columns to FEATURE_1..N
    merged = _renumber_features(merged)

    saved_path: Union[str, None] = None
    if output_path is not None:
        merged.to_csv(output_path, index=False)
        saved_path = str(output_path)

    return merged, saved_path
