#!/usr/bin/env python3
"""Compatibility wrapper for active SPECTRAPRINTS_NMR_V2 preparation."""

from __future__ import annotations

import csv
import os
from pathlib import Path

from .preparation import (
    PreparationTask,
    _worker_init,
    prepare_batch,
    prepare_mol_v2,
    prepare_task,
)


def generate_mol_files(csv_path: str, strict_mode: bool = True) -> str:
    """Generate deterministic NMR V2 MOL files for the historical API.

    ``strict_mode`` is accepted only for API compatibility. NMR V2 has one
    frozen preparation path and no ETKDG/coordinate-origin branch.
    """
    del strict_mode
    output = Path.cwd() / "mols"
    output.mkdir(parents=True, exist_ok=True)
    records = []
    with Path(csv_path).open(newline="", encoding="utf-8") as handle:
        for index, row in enumerate(csv.DictReader(handle)):
            name = str(row.get("MOLECULE_NAME") or f"m{index:08d}")
            records.append({
                "internal_id": name,
                "molecule_name": name,
                "smiles": str(row["SMILES"]),
            })
    results = prepare_batch(records, output)
    failures = [result for result in results if not result.successful]
    if failures:
        error_log = Path.cwd() / "mol_creation_error.log"
        error_log.write_text(
            "==== NMR V2 MOL CREATION ERRORS ====\n\n" + "\n".join(
                f"Molecule: {item.molecule_name}\nSMILES: {item.smiles}\n"
                f"Error: {item.error_type}: {item.error_message}\n"
                for item in failures
            ),
            encoding="utf-8",
        )
    return os.fspath(output)


__all__ = [
    "PreparationTask",
    "_worker_init",
    "prepare_batch",
    "prepare_mol_v2",
    "prepare_task",
    "generate_mol_files",
]
