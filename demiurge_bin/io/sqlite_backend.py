"""Streaming SQLite input and compact float32-vector output."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from demiurge_bin.run_state import atomic_write_json, atomic_write_text, file_sha256, utc_now

from .base import ColumnSelector, InputDescription, InputReader, OutputWriter, record_from_values, resolve_column
from .csv_backend import _aggregate_failures, _retain_artifacts, _validate_replay


SQLITE_SCHEMA_VERSION = 1
FEATURE_DTYPE = "<f4"


def quote_identifier(value: str) -> str:
    name = str(value).strip()
    if not name or "\x00" in name:
        raise ValueError("SQLite identifier must be non-empty and contain no NUL")
    return '"' + name.replace('"', '""') + '"'


def _connect_readonly(path: Path) -> sqlite3.Connection:
    uri = f"file:{path.as_posix()}?mode=ro&immutable=1"
    connection = sqlite3.connect(uri, uri=True)
    connection.execute("PRAGMA query_only=ON")
    return connection


def _reject_live_sidecars(path: Path) -> None:
    for suffix in ("-wal", "-shm", "-journal"):
        sidecar = path.with_name(path.name + suffix)
        if sidecar.is_file() and sidecar.stat().st_size:
            raise RuntimeError(
                f"SQLite input has an active {suffix} sidecar; checkpoint the database before use: {sidecar}"
            )


class SqliteInputReader(InputReader):
    def __init__(
        self,
        path: Path,
        *,
        table: str | None,
        query: str | None,
        id_column: str,
        smiles_column: str,
        label_column: ColumnSelector,
    ):
        if bool(table) == bool(query):
            raise ValueError("SQLite input requires exactly one of input_table or input_query")
        self.path = path
        self.table = table
        self.query = self._validated_query(query) if query else None
        self.id_column = id_column
        self.smiles_column = smiles_column
        self.label_column = label_column
        self._description: InputDescription | None = None
        self._label_name: str | None = None
        self._select_sql: str | None = None

    @staticmethod
    def _validated_query(query: str) -> str:
        selected = query.strip()
        if selected.endswith(";"):
            selected = selected[:-1].rstrip()
        if not selected or not selected.lstrip().lower().startswith(("select ", "with ")):
            raise ValueError("input_query must be one read-only SELECT statement")
        if ";" in selected:
            raise ValueError("input_query must contain exactly one statement")
        return selected

    def _base_sql(self) -> str:
        if self.query is not None:
            return f"({self.query})"
        assert self.table is not None
        return quote_identifier(self.table)

    def describe(self) -> InputDescription:
        if self._description is not None:
            return self._description
        _reject_live_sidecars(self.path)
        with closing(_connect_readonly(self.path)) as connection:
            if self.table is not None:
                exists = connection.execute(
                    "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?",
                    (self.table,),
                ).fetchone()
                if exists is None:
                    raise ValueError(f"SQLite input table/view does not exist: {self.table!r}")
                probe = connection.execute(f"SELECT * FROM {quote_identifier(self.table)} LIMIT 0")
            else:
                probe = connection.execute(f"SELECT * FROM ({self.query}) LIMIT 0")
            columns = [str(item[0]) for item in probe.description or ()]
            for option, name in (("id_column", self.id_column), ("smiles_column", self.smiles_column)):
                if columns.count(name) != 1:
                    reason = "does not exist" if name not in columns else "is ambiguous"
                    raise ValueError(f"{option}={name!r} {reason} in the SQLite input")
            _, self._label_name = resolve_column(columns, self.label_column, "label_column")
            selected = ", ".join(quote_identifier(name) for name in (
                self.id_column, self.smiles_column, self._label_name,
            ))
            self._select_sql = (
                f"SELECT {selected} FROM {self._base_sql()} "
                f"ORDER BY {quote_identifier(self.id_column)}"
            )
            total = int(connection.execute(f"SELECT COUNT(*) FROM {self._base_sql()}").fetchone()[0])
        self._description = InputDescription(
            total_rows=total,
            label_name=self._label_name,
            configuration={
                "input_format": "sqlite",
                "input_table": self.table,
                "input_query": self.query,
                "id_column": self.id_column,
                "smiles_column": self.smiles_column,
                "label_column": self.label_column,
                "label_name": self._label_name,
                "ordering": self.id_column,
            },
        )
        return self._description

    def iter_batches(self, batch_size: int, start_index: int = 0) -> Iterator[tuple[int, list[dict[str, Any]]]]:
        self.describe()
        assert self._select_sql is not None and self._label_name is not None
        _reject_live_sidecars(self.path)
        with closing(_connect_readonly(self.path)) as connection:
            cursor = connection.execute(self._select_sql + " LIMIT -1 OFFSET ?", (int(start_index),))
            source_index = int(start_index)
            while True:
                rows = cursor.fetchmany(batch_size)
                if not rows:
                    return
                offset = source_index
                records = []
                for molecule_id, smiles, label in rows:
                    records.append(record_from_values(
                        source_index, molecule_id, smiles, label, self._label_name, self.id_column
                    ))
                    source_index += 1
                yield offset, records


class SqliteOutputWriter(OutputWriter):
    def __init__(
        self,
        output_root: Path,
        output_db: Path,
        *,
        output_table: str,
        metadata_table: str,
        feature_contract: dict[str, Any],
        overwrite: bool,
    ):
        self.output_root = output_root
        self.output_db = output_db
        self.output_table = output_table
        self.metadata_table = metadata_table
        self.feature_contract = feature_contract
        self.overwrite = overwrite
        quote_identifier(output_table)
        quote_identifier(metadata_table)
        if output_table == metadata_table:
            raise ValueError("output_table and metadata_table must differ")

    def configuration(self) -> dict[str, Any]:
        return {
            "output_format": "sqlite",
            "output_db": str(self.output_db),
            "output_table": self.output_table,
            "metadata_table": self.metadata_table,
            "feature_blob_dtype": FEATURE_DTYPE,
        }

    def initialize(self, *, resume: bool) -> None:
        self.output_db.parent.mkdir(parents=True, exist_ok=True)
        if resume:
            if not self.output_db.is_file():
                raise RuntimeError(f"SQLite resume output is missing: {self.output_db}")
            with closing(sqlite3.connect(self.output_db)) as connection:
                for table in (self.output_table, self.metadata_table):
                    exists = connection.execute(
                        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
                    ).fetchone()
                    if exists is None:
                        raise RuntimeError(f"SQLite resume table is missing: {table!r}")
            return
        if self.output_db.exists():
            if not self.overwrite:
                raise FileExistsError(
                    f"SQLite output already exists; choose a new path or pass --overwrite-output: {self.output_db}"
                )
            self.output_db.unlink()
        with closing(sqlite3.connect(self.output_db)) as connection:
            connection.execute("PRAGMA journal_mode=DELETE")
            connection.execute("PRAGMA synchronous=FULL")
            connection.execute(f"""
                CREATE TABLE {quote_identifier(self.output_table)} (
                    source_index INTEGER PRIMARY KEY,
                    molecule_id TEXT NOT NULL,
                    label,
                    status TEXT NOT NULL CHECK(status IN ('SUCCESS','FAILED')),
                    error_stage TEXT,
                    error_type TEXT,
                    error TEXT,
                    feature_blob BLOB,
                    feature_sha256 TEXT,
                    CHECK((status='SUCCESS' AND feature_blob IS NOT NULL) OR
                          (status='FAILED' AND feature_blob IS NULL))
                )
            """)
            connection.execute(f"""
                CREATE TABLE {quote_identifier(self.metadata_table)} (
                    key TEXT PRIMARY KEY,
                    value_json TEXT NOT NULL
                )
            """)
            contract_id = self.feature_contract["contract_id"]
            component_counts = {
                "DEMIURGE_1H_NMR_V2": (200, 0, 0),
                "DEMIURGE_13C_NMR_V2": (0, 200, 0),
                "DEMIURGE_HYBRID_NMR_V2_H_C": (200, 200, 0),
                "DEMIURGE_ECFP4": (0, 0, 2048),
                "DEMIURGE_TOTAL_NMR_V2_H_C_ECFP4": (200, 200, 2048),
            }
            h_count, c_count, fp_count = component_counts[contract_id]
            metadata = {
                "schema_version": SQLITE_SCHEMA_VERSION,
                "scientific_core_version": contract_id,
                "feature_representation": contract_id,
                "nmr_representation_version": self.feature_contract.get("nmr_representation_version"),
                "feature_count": int(self.feature_contract["feature_dimension"]),
                "feature_blob_bytes": int(self.feature_contract["feature_dimension"]) * np.dtype(FEATURE_DTYPE).itemsize,
                "dtype": FEATURE_DTYPE,
                "byte_order": "little-endian",
                "1h_feature_count": h_count,
                "13c_feature_count": c_count,
                "ecfp4_feature_count": fp_count,
                "feature_offsets": {
                    "1H": [0, h_count] if h_count else None,
                    "13C": [h_count, h_count + c_count] if c_count else None,
                    "ECFP4": [h_count + c_count, h_count + c_count + fp_count] if fp_count else None,
                },
                "feature_order": self.feature_contract["feature_order"],
                "created_at": utc_now(),
            }
            connection.executemany(
                f"INSERT INTO {quote_identifier(self.metadata_table)}(key,value_json) VALUES (?,?)",
                [(key, json.dumps(value, sort_keys=True, ensure_ascii=False)) for key, value in metadata.items()],
            )
            connection.commit()

    def _encoded_rows(self, feature_rows: list[dict[str, Any]], metadata: list[dict[str, Any]]) -> list[tuple[Any, ...]]:
        features = {int(item["source_index"]): item for item in feature_rows}
        encoded = []
        expected_bytes = int(self.feature_contract["feature_dimension"]) * np.dtype(FEATURE_DTYPE).itemsize
        for item in sorted(metadata, key=lambda value: int(value["source_index"])):
            source_index = int(item["source_index"])
            feature = features.get(source_index)
            if item["status"] == "SUCCESS":
                if feature is None:
                    raise RuntimeError(f"Successful metadata has no feature row: source_index={source_index}")
                blob = np.asarray(feature["vector"], dtype=FEATURE_DTYPE).tobytes(order="C")
                if len(blob) != expected_bytes:
                    raise RuntimeError(f"Unexpected feature blob length at source_index={source_index}")
                digest = hashlib.sha256(blob).hexdigest()
            else:
                blob = None
                digest = None
            encoded.append((
                source_index,
                item["molecule_name"],
                feature["label"] if feature is not None else item.get("label"),
                item["status"],
                item.get("failure_stage"),
                item.get("failure_type"),
                item.get("failure_message"),
                blob,
                digest,
            ))
        return encoded

    def commit_batch(
        self,
        batch_index: int,
        feature_rows: list[dict[str, Any]],
        metadata: list[dict[str, Any]],
        timings: dict[str, float],
        scratch_batch: Path,
        retain_scientific_artifacts: bool,
    ) -> Path:
        encoded = self._encoded_rows(feature_rows, metadata)
        table = quote_identifier(self.output_table)
        with closing(sqlite3.connect(self.output_db)) as connection:
            connection.execute("PRAGMA synchronous=FULL")
            for row in encoded:
                existing = connection.execute(
                    f"SELECT source_index,molecule_id,label,status,error_stage,error_type,error,feature_blob,feature_sha256 FROM {table} WHERE source_index=?",
                    (row[0],),
                ).fetchone()
                if existing is not None:
                    if existing != row:
                        raise RuntimeError(f"Existing SQLite row differs during replay: source_index={row[0]}")
                    continue
                connection.execute(
                    f"INSERT INTO {table}(source_index,molecule_id,label,status,error_stage,error_type,error,feature_blob,feature_sha256) VALUES (?,?,?,?,?,?,?,?,?)",
                    row,
                )
            connection.commit()
        batches = self.output_root / "batches"
        batches.mkdir(parents=True, exist_ok=True)
        final = batches / f"batch_{batch_index:08d}"
        temporary = batches / f".batch_{batch_index:08d}.{os.getpid()}.tmp"
        if temporary.exists():
            shutil.rmtree(temporary)
        temporary.mkdir()
        metadata_path = temporary / "metadata.jsonl"
        atomic_write_text(metadata_path, "".join(
            json.dumps(item, sort_keys=True, ensure_ascii=False) + "\n" for item in metadata
        ))
        _retain_artifacts(temporary, scratch_batch, retain_scientific_artifacts)
        atomic_write_json(temporary / "commit.json", {
            "schema_version": 1,
            "batch_index": batch_index,
            "successful": len(feature_rows),
            "failed": sum(1 for item in metadata if item["status"] == "FAILED"),
            "metadata_sha256": file_sha256(metadata_path),
            "sqlite_output": str(self.output_db),
            "timing_seconds": timings,
            "committed_at": utc_now(),
        })
        if final.exists():
            _validate_replay(final, temporary, retain_scientific_artifacts, include_features=False)
            shutil.rmtree(temporary)
            return final
        os.replace(temporary, final)
        return final

    def finalize(self, *, total: int, successful: int, failed: int) -> tuple[Path, Path]:
        with closing(sqlite3.connect(self.output_db)) as connection:
            table = quote_identifier(self.metadata_table)
            values = {
                "total_rows": total,
                "successful_rows": successful,
                "failed_rows": failed,
                "completed_at": utc_now(),
            }
            connection.executemany(
                f"INSERT OR REPLACE INTO {table}(key,value_json) VALUES (?,?)",
                [(key, json.dumps(value, sort_keys=True)) for key, value in values.items()],
            )
            connection.commit()
        return self.output_db, _aggregate_failures(self.output_root)
