"""Authoritative SPECTRAPRINTS_NMR_V2 molecule preparation for Demiurge."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


def require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:  # pragma: no cover - dependency preflight
        raise RuntimeError("NMR V2 molecule preparation requires RDKit") from exc
    return Chem, AllChem


def prepare_mol_v2(smiles: str) -> tuple[str, str]:
    """Return canonical SMILES and deterministic explicit-H V3000 2D MOL."""
    Chem, AllChem = require_rdkit()
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError("RDKit rejected SMILES")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    prepared = Chem.AddHs(Chem.MolFromSmiles(canonical))
    if prepared is None:
        raise ValueError("RDKit failed to recreate canonical SMILES")
    AllChem.Compute2DCoords(prepared)
    Chem.RemoveStereochemistry(prepared)
    block = Chem.MolToMolBlock(prepared, forceV3000=True)
    if "V3000" not in block:
        raise RuntimeError("RDKit did not produce V3000 output")
    return canonical, block


def _worker_init() -> None:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")


@dataclass(frozen=True)
class PreparationTask:
    internal_id: str
    molecule_name: str
    smiles: str
    output_directory: str


@dataclass(frozen=True)
class PreparationResult:
    internal_id: str
    molecule_name: str
    smiles: str
    canonical_smiles: str | None
    mol_path: str | None
    error_type: str | None
    error_message: str | None

    @property
    def successful(self) -> bool:
        return self.mol_path is not None


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def prepare_task(task: PreparationTask) -> PreparationResult:
    try:
        canonical, block = prepare_mol_v2(task.smiles)
        output = Path(task.output_directory) / f"{task.internal_id}.mol"
        _atomic_text(output, block)
        return PreparationResult(
            task.internal_id,
            task.molecule_name,
            task.smiles,
            canonical,
            str(output),
            None,
            None,
        )
    except Exception as exc:  # molecule-level scientific/QC result
        return PreparationResult(
            task.internal_id,
            task.molecule_name,
            task.smiles,
            None,
            None,
            type(exc).__name__,
            str(exc),
        )


def prepare_batch(
    records: Iterable[dict[str, Any]],
    output_directory: Path,
    *,
    pool: Any = None,
) -> list[PreparationResult]:
    output_directory.mkdir(parents=True, exist_ok=True)
    tasks = [
        PreparationTask(
            str(record["internal_id"]),
            str(record["molecule_name"]),
            str(record["smiles"]),
            str(output_directory),
        )
        for record in records
    ]
    results = list(map(prepare_task, tasks)) if pool is None else list(pool.map(prepare_task, tasks))
    return sorted(results, key=lambda result: result.internal_id)
