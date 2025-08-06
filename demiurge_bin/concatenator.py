#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
concatenator.py
===============

Utility for merging two machine-learning input tables
(¹H and ¹³C pseudo-spectra) into a single *hybrid* dataset.

Both inputs must contain two metadata columns:

    • MOLECULE_NAME
    • LABEL                (target value)

All remaining columns are treated as numeric features and will be
prefixed with ``H_`` or ``C_`` before concatenation.  The final table
layout is::

    MOLECULE_NAME, LABEL, FEATURE_1, FEATURE_2, …

Example
-------

    from concatenator import concatenate

    hybrid_df, _ = concatenate(
        ["spectra_1H.csv", "spectra_13C.csv"],
        output_path="hybrid.csv",
    )
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple, Union

import pandas as pd

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
META_COLS: list[str] = ["MOLECULE_NAME", "LABEL"]

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def _to_dataframe(data: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    """Return *data* as a fresh ``pandas.DataFrame``."""
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return pd.read_csv(data)


def _split_meta_features(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split *df* into metadata and feature frames."""
    meta_df = df[META_COLS]
    feat_df = df.drop(columns=META_COLS)
    return meta_df, feat_df


def _renumber_features(df: pd.DataFrame) -> pd.DataFrame:
    """Rename every feature column to ``FEATURE_n`` while preserving order."""
    feature_cols = [c for c in df.columns if c not in META_COLS]
    mapping = {old: f"FEATURE_{i + 1}" for i, old in enumerate(feature_cols)}
    return df.rename(columns=mapping)

# -----------------------------------------------------------------------------
# Core
# -----------------------------------------------------------------------------
def _concat_features(df_h: pd.DataFrame, df_c: pd.DataFrame) -> pd.DataFrame:
    """Merge ¹H and ¹³C feature matrices on *META_COLS* and renumber columns."""
    meta_h, feats_h = _split_meta_features(df_h)
    meta_c, feats_c = _split_meta_features(df_c)

    feats_h.columns = [f"H_{col}" for col in feats_h.columns]
    feats_c.columns = [f"C_{col}" for col in feats_c.columns]

    df_h_full = pd.concat([meta_h, feats_h], axis=1)
    df_c_full = pd.concat([meta_c, feats_c], axis=1)

    merged = pd.merge(df_h_full, df_c_full, on=META_COLS, how="inner")
    return _renumber_features(merged)


def concatenate(
    datasets: Iterable[Union[str, Path, pd.DataFrame]],
    output_path: Union[str, Path, None] = None,
) -> Tuple[pd.DataFrame, Union[str, None]]:
    """
    Concatenate two spectra tables (**H** and **C**) into a hybrid table.

    Parameters
    ----------
    datasets
        Iterable with exactly two elements - each either a CSV path or a
        ``pandas.DataFrame``.
    output_path
        Optional path for saving the merged CSV.  If *None*, the file is not
        written.

    Returns
    -------
    (merged_df, saved_path)
        merged_df : ``pandas.DataFrame``
            The combined dataset.
        saved_path : str | None
            Location where the CSV was saved, or *None* if nothing was saved.
    """
    ds_list = list(datasets)
    if len(ds_list) != 2:
        raise ValueError("Expected exactly two datasets: one ¹H and one ¹³C.")

    df_h, df_c = map(_to_dataframe, ds_list)
    merged_df = _concat_features(df_h, df_c)

    saved_path: Union[str, None] = None
    if output_path is not None:
        merged_df.to_csv(output_path, index=False)
        saved_path = str(output_path)

    return merged_df, saved_path
