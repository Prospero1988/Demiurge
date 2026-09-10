"""Construct validated I/O backends without exposing formats to science code."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .base import ColumnSelector, InputReader, OutputWriter
from .csv_backend import CsvInputReader, CsvOutputWriter
from .sqlite_backend import SqliteInputReader, SqliteOutputWriter


def create_input_reader(
    path: Path,
    *,
    input_format: str,
    input_table: str | None,
    input_query: str | None,
    id_column: str,
    smiles_column: str,
    label_column: ColumnSelector,
) -> InputReader:
    if input_format == "csv":
        if input_table or input_query:
            raise ValueError("input_table/input_query are valid only for SQLite input")
        return CsvInputReader(path, id_column, smiles_column, label_column)
    if input_format == "sqlite":
        return SqliteInputReader(
            path,
            table=input_table,
            query=input_query,
            id_column=id_column,
            smiles_column=smiles_column,
            label_column=label_column,
        )
    raise ValueError(f"Unsupported input_format: {input_format!r}")


def create_output_writer(
    *,
    output_format: str,
    output_root: Path,
    output_db: Path | None,
    output_table: str,
    metadata_table: str,
    input_stem: str,
    mode: str,
    feature_contract: dict[str, Any],
    overwrite_output: bool,
) -> OutputWriter:
    if output_format == "csv":
        if output_db is not None:
            raise ValueError("output_db is valid only for SQLite output")
        if overwrite_output:
            raise ValueError("overwrite_output is valid only for SQLite output")
        if output_table != "demiurge_features" or metadata_table != "demiurge_metadata":
            raise ValueError("output_table/metadata_table are valid only for SQLite output")
        return CsvOutputWriter(output_root, input_stem, mode, int(feature_contract["feature_dimension"]))
    if output_format == "sqlite":
        selected_db = output_db or (
            output_root / "generated_ML_inputs" / f"{input_stem}_{mode}_ML_input.sqlite"
        )
        return SqliteOutputWriter(
            output_root,
            selected_db,
            output_table=output_table,
            metadata_table=metadata_table,
            feature_contract=feature_contract,
            overwrite=overwrite_output,
        )
    raise ValueError(f"Unsupported output_format: {output_format!r}")
