"""Scratch-backed molecular feature cache for record-preserving de-duplication."""

from __future__ import annotations

import json
import shutil
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class CachedFeatures:
    identity_smiles: str
    vector: list[int]
    diagnostics: dict[str, Any]


class ComputationCache:
    """Store one vector per molecular identity without retaining it in RAM."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._artifacts = path.parent / "molecular_artifact_cache"
        self._artifacts.mkdir(parents=True, exist_ok=True)
        self._connection = sqlite3.connect(path)
        self._connection.execute("PRAGMA journal_mode=WAL")
        self._connection.execute("PRAGMA synchronous=FULL")
        self._connection.execute(
            """CREATE TABLE IF NOT EXISTS feature_cache (
                identity_sha256 TEXT PRIMARY KEY,
                identity_smiles TEXT NOT NULL,
                vector_int32 BLOB NOT NULL,
                vector_length INTEGER NOT NULL,
                diagnostics_json TEXT NOT NULL
            )"""
        )
        self._connection.execute(
            """CREATE TABLE IF NOT EXISTS accepted_records (
                record_id TEXT PRIMARY KEY,
                molecule_name TEXT NOT NULL,
                raw_smiles TEXT NOT NULL,
                label_json TEXT NOT NULL,
                identity_sha256 TEXT NOT NULL
            )"""
        )
        self._connection.commit()

    def get(self, identity_sha256: str, identity_smiles: str) -> CachedFeatures | None:
        row = self._connection.execute(
            "SELECT identity_smiles,vector_int32,vector_length,diagnostics_json "
            "FROM feature_cache WHERE identity_sha256=?",
            (identity_sha256,),
        ).fetchone()
        if row is None:
            return None
        stored_smiles, blob, length, diagnostics = row
        if stored_smiles != identity_smiles:
            raise RuntimeError("Molecular identity SHA256 collision detected")
        vector = np.frombuffer(blob, dtype="<i4")
        if len(vector) != int(length):
            raise RuntimeError("Corrupt molecular feature cache entry")
        return CachedFeatures(stored_smiles, vector.astype(int).tolist(), json.loads(diagnostics))

    def put(
        self,
        identity_sha256: str,
        identity_smiles: str,
        vector: list[int],
        diagnostics: dict[str, Any],
    ) -> None:
        blob = np.asarray(vector, dtype="<i4").tobytes(order="C")
        encoded = json.dumps(diagnostics, sort_keys=True, separators=(",", ":"))
        prior = self._connection.execute(
            "SELECT identity_smiles,vector_int32,vector_length,diagnostics_json "
            "FROM feature_cache WHERE identity_sha256=?",
            (identity_sha256,),
        ).fetchone()
        candidate = (identity_smiles, blob, len(vector), encoded)
        if prior is not None:
            if prior != candidate:
                raise RuntimeError("Molecular feature cache replay differs")
            return
        self._connection.execute(
            "INSERT INTO feature_cache VALUES (?,?,?,?,?)",
            (identity_sha256, identity_smiles, blob, len(vector), encoded),
        )
        self._connection.commit()

    def close(self) -> None:
        self._connection.close()

    def capture_artifacts(
        self,
        identity_sha256: str,
        *,
        mol_path: Path | None,
        raw_h_path: Path | None,
        raw_c_path: Path | None,
    ) -> None:
        destination = self._artifacts / identity_sha256
        destination.mkdir(parents=True, exist_ok=True)
        for name, source in (("molecule.mol", mol_path), ("1H.csv", raw_h_path), ("13C.csv", raw_c_path)):
            if source is not None and source.is_file():
                target = destination / name
                if target.exists() and target.read_bytes() != source.read_bytes():
                    raise RuntimeError("Molecular artifact cache replay differs")
                if not target.exists():
                    shutil.copyfile(source, target)

    def restore_artifacts(
        self,
        identity_sha256: str,
        internal_id: str,
        *,
        mol_directory: Path,
        raw_h_directory: Path,
        raw_c_directory: Path,
    ) -> None:
        source = self._artifacts / identity_sha256
        for name, target in (
            ("molecule.mol", mol_directory / f"{internal_id}.mol"),
            ("1H.csv", raw_h_directory / f"{internal_id}.csv"),
            ("13C.csv", raw_c_directory / f"{internal_id}.csv"),
        ):
            cached = source / name
            if cached.is_file():
                shutil.copyfile(cached, target)

    def observe_record(
        self,
        *,
        record_id: str,
        molecule_name: str,
        raw_smiles: str,
        label: Any,
        identity_sha256: str,
    ) -> None:
        candidate = (
            record_id, molecule_name, raw_smiles,
            json.dumps(label, sort_keys=True, ensure_ascii=False), identity_sha256,
        )
        prior = self._connection.execute(
            "SELECT record_id,molecule_name,raw_smiles,label_json,identity_sha256 "
            "FROM accepted_records WHERE record_id=?", (record_id,),
        ).fetchone()
        if prior is not None:
            if prior != candidate:
                raise RuntimeError("Accepted-record audit replay differs")
            return
        self._connection.execute("INSERT INTO accepted_records VALUES (?,?,?,?,?)", candidate)
        self._connection.commit()

    def audit(self) -> dict[str, int]:
        total, names, identities = self._connection.execute(
            "SELECT COUNT(*),COUNT(DISTINCT molecule_name),COUNT(DISTINCT identity_sha256) "
            "FROM accepted_records"
        ).fetchone()
        exact_duplicates = self._connection.execute(
            "SELECT COALESCE(SUM(count_rows - 1),0) FROM ("
            "SELECT COUNT(*) AS count_rows FROM accepted_records "
            "GROUP BY molecule_name,raw_smiles,label_json HAVING COUNT(*) > 1)"
        ).fetchone()[0]
        computations = self._connection.execute("SELECT COUNT(*) FROM feature_cache").fetchone()[0]
        return {
            "accepted_output_records": int(total),
            "unique_molecule_names": int(names),
            "duplicate_molecule_name_records": int(total - names),
            "unique_molecular_structures": int(identities),
            "reused_molecular_structure_records": int(total - identities),
            "exact_duplicate_input_records": int(exact_duplicates),
            "unique_feature_computations": int(computations),
        }
