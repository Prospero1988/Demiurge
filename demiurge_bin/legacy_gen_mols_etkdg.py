# ARCHIVED REFERENCE ONLY. Not imported by the production NMR V2 pipeline.
#!/usr/bin/env python3
"""
Generate 2D MOL files from SMILES strings (parallel, chunked by CPU).

Pipeline
--------
1.  Canonicalise the SMILES.
2.  Add explicit hydrogens.
3.  Try to embed a 3D conformer with ETKDG (3 retries).
    • If ETKDG fails → fall back to RDKit CoordGen (2D).
4.  Force a switch to OpenBabel for molecules that
    contain hyper-valent sulphur (valence > 4) or after an
    ETKDG failure.
    • obabel -d --gen2D strips wedge bonds and flattens the structure.
5.  Write a V3000 MOL file (2D coordinates, no stereo wedges).
6.  Log errors to *mol_creation_error.log* and all fall-backs/
    warnings to *mol_creation_warning.log*.

Notes
-----
*   OpenBabel (`obabel`) must be on your system `PATH`
    (e.g. `conda install -c conda-forge openbabel`).
*   RDKit warnings are silenced via `RDLogger.DisableLog('rdApp.*')`.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from multiprocessing import Pool
from typing import List, Tuple, Optional

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdCoordGen


# ──────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────
ANSI_GREEN = "\033[38;5;46m"
ANSI_RED = "\033[38;5;196m"
ANSI_ORANGE = "\033[38;5;214m"
ANSI_RESET = "\033[0m"

PROGRESS_BAR_LEN = 25
MAX_ETKDG_RETRIES = 3
EMBED_RANDOM_SEED = 42


# ──────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────
def _worker_init() -> None:
    """Initializer for worker processes."""
    RDLogger.DisableLog("rdApp.*")


def safe_embed_molecule(
    mol: Chem.Mol,
    max_retries: int = MAX_ETKDG_RETRIES,
    seed: int = EMBED_RANDOM_SEED,
) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """
    Try ETKDG embedding up to *max_retries* times, fall back to CoordGen.

    Returns
    -------
    mol
        RDKit molecule with at least one conformer, or None on hard fail.
    warning
        None on ETKDG success, otherwise a human-readable note.
    """
    try:
        params = AllChem.ETKDGv3()
    except AttributeError:
        try:
            params = AllChem.ETKDGv2()
        except AttributeError:
            params = AllChem.ETKDG()

    params.randomSeed = seed

    for _ in range(max_retries):
        mol.RemoveAllConformers()
        if AllChem.EmbedMolecule(mol, params) == 0:
            return mol, None

    try:
        rdCoordGen.AddCoords(mol)
        return mol, f"ETKDG failed ({max_retries}x) → used CoordGen"
    except Exception as exc:
        return None, f"ETKDG + CoordGen failed: {exc}"


def canonical_smiles(smiles: str) -> str:
    """Return RDKit-canonical SMILES or raise ValueError."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


def openbabel_fallback(rdkit_mol: Chem.Mol, out_path: str) -> Tuple[bool, Optional[str]]:
    """
    Run obabel -d --gen2D on rdkit_mol; write to out_path.
    Returns (success, error_message).
    """
    with tempfile.NamedTemporaryFile(suffix=".mol", delete=False) as tmp:
        tmp.write(Chem.MolToMolBlock(rdkit_mol, forceV3000=True).encode())
        tmp_path = tmp.name

    cmd = ["obabel", tmp_path, "-O", out_path, "-d", "--gen2D"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.remove(tmp_path)
        return True, None
    except subprocess.CalledProcessError as exc:
        os.remove(tmp_path)
        return False, exc.stderr.decode().strip()


def print_progress(current: int, total: int) -> None:
    """Draw a coloured, in-place ASCII progress bar."""
    if total <= 0:
        return
    filled = int(PROGRESS_BAR_LEN * current / total)
    bar = ANSI_GREEN + "█" * filled + "-" * (PROGRESS_BAR_LEN - filled) + ANSI_RESET
    percent = int(100 * current / total)
    sys.stdout.write(f"\rProgress: |{bar}| {current}/{total} ({percent}%)")
    sys.stdout.flush()
    if current >= total:
        print()


# ──────────────────────────────────────────────────────────────
# Extra “weird-chemistry” detectors
# ──────────────────────────────────────────────────────────────
EXOTIC_VALENCE_LIMITS = {"S": 4, "P": 4, "As": 4, "Se": 4}
TRANSITION_METALS = {
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    57, 72, 73, 74, 75, 76, 77, 78, 79
}


def needs_openbabel(mol: Chem.Mol) -> bool:
    """
    Return True for molecules that should bypass RDKit-only handling.
    Triggers:
        • hyper-valent S/P/As/Se
        • any transition metal
        • radicals
        • too many heavy atoms (> 150)
    """
    if any(
        a.GetSymbol() in EXOTIC_VALENCE_LIMITS
        and a.GetTotalValence() > EXOTIC_VALENCE_LIMITS[a.GetSymbol()]
        for a in mol.GetAtoms()
    ):
        return True

    if any(a.GetAtomicNum() in TRANSITION_METALS for a in mol.GetAtoms()):
        return True

    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return True

    if mol.GetNumHeavyAtoms() > 150:
        return True

    return False


def is_dot_smiles(smiles: str) -> bool:
    """True if SMILES contains disconnected fragments (“dot-SMILES”)."""
    return "." in smiles


# ──────────────────────────────────────────────────────────────
# Parallel processing: chunk worker
# ──────────────────────────────────────────────────────────────
@dataclass
class ChunkResult:
    processed: int
    saved: int
    errors: List[str]
    warnings: List[str]


def process_chunk(args: Tuple[List[Tuple[str, str]], str, bool]) -> ChunkResult:
    """
    Process a list of (name, raw_smiles) tuples and write MOL files.
    Returns aggregated logs.
    """
    rows, output_dir, strict_mode = args
    errors: List[str] = []
    warnings: List[str] = []
    saved_files = 0

    for name, raw_smiles in rows:
        try:
            smiles = canonical_smiles(raw_smiles)
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))

            mol, warn_msg = safe_embed_molecule(mol)
            if mol is None:
                raise ValueError(warn_msg)
            if warn_msg:
                warnings.append(f"{name}: {warn_msg}")

            force_babel = False
            reasons: List[str] = []

            if needs_openbabel(mol):
                force_babel = True
                reasons.append("exotic atom / metal / radical / size")

            if is_dot_smiles(smiles):
                force_babel = True
                reasons.append("dot-SMILES (disconnected fragments)")

            if warn_msg:
                force_babel = True
                reasons.append("ETKDG failure")

            if force_babel:
                warnings.append(f"{name}: OpenBabel fallback → {', '.join(reasons)}")

            conf = mol.GetConformer()
            if strict_mode and all(
                conf.GetAtomPosition(i).Length() < 0.1 for i in range(mol.GetNumAtoms())
            ):
                raise ValueError("All atoms at origin (invalid 3D)")

            mol2d = Chem.Mol(mol)
            AllChem.Compute2DCoords(mol2d)
            Chem.RemoveStereochemistry(mol2d)

            out_path = os.path.join(output_dir, f"{name}.mol")

            if not force_babel:
                with open(out_path, "w", encoding="utf-8") as handle:
                    handle.write(Chem.MolToMolBlock(mol2d, forceV3000=True))
            else:
                success, ob_error = openbabel_fallback(mol2d, out_path)
                if not success:
                    raise ValueError(f"OpenBabel fallback failed: {ob_error}")

            saved_files += 1

        except Exception as exc:  # pylint: disable=broad-except
            errors.append(
                f"Molecule: {name}\nSMILES: {raw_smiles}\nError: {exc}\n"
            )

    return ChunkResult(processed=len(rows), saved=saved_files, errors=errors, warnings=warnings)


def split_into_n_chunks(items: List[Tuple[str, str]], n: int) -> List[List[Tuple[str, str]]]:
    """Split items into n chunks as evenly as possible."""
    if n <= 1:
        return [items]
    k, m = divmod(len(items), n)
    chunks: List[List[Tuple[str, str]]] = []
    start = 0
    for i in range(n):
        size = k + (1 if i < m else 0)
        chunks.append(items[start:start + size])
        start += size
    return [c for c in chunks if c]  # drop empty chunks


# ──────────────────────────────────────────────────────────────
# Main routine
# ──────────────────────────────────────────────────────────────
def generate_mol_files(csv_path: str, strict_mode: bool = True) -> str:
    """
    Convert SMILES in csv_path to flat MOL files (parallel, CPU-chunked).

    Parameters
    ----------
    csv_path
        CSV with columns MOLECULE_NAME and SMILES.
    strict_mode
        If True, reject molecules whose coords all sit near (0, 0, 0).

    Returns
    -------
    str
        Output directory path.
    """
    output_dir = os.path.join(os.getcwd(), "mols")
    os.makedirs(output_dir, exist_ok=True)

    data = pd.read_csv(csv_path)
    data = data.drop_duplicates(subset="MOLECULE_NAME", keep="first")

    rows: List[Tuple[str, str]] = [(r.MOLECULE_NAME, r.SMILES) for r in data.itertuples(index=False)]
    total = len(rows)

    if total == 0:
        print(f"{ANSI_ORANGE}No rows found in CSV. Nothing to do.{ANSI_RESET}")
        return output_dir

    cpu_count = os.cpu_count() or 1
    n_workers = min(cpu_count, total)

    chunk_factor = 4
    n_chunks = min(total, n_workers * chunk_factor)
    chunks = split_into_n_chunks(rows, n_chunks)

    print("\nGenerating *.mol files …\n")
    print(f"Detected CPUs: {cpu_count} → using workers: {n_workers} (chunks: {len(chunks)})\n")

    print_progress(0, total)

    all_errors: List[str] = []
    all_warnings: List[str] = []
    saved_files = 0
    done = 0

    # Prepare args for each worker: (chunk_rows, output_dir, strict_mode)
    worker_args = [(chunk, output_dir, strict_mode) for chunk in chunks]

    with Pool(processes=n_workers, initializer=_worker_init) as pool:
        for result in pool.imap_unordered(process_chunk, worker_args):
            done += result.processed
            saved_files += result.saved
            all_errors.extend(result.errors)
            all_warnings.extend(result.warnings)
            print_progress(done, total)

    # ── Write logs ─────────────────────────────────────────────────────
    if all_errors:
        with open("mol_creation_error.log", "w", encoding="utf-8") as fh_err:
            fh_err.write("==== MOL CREATION ERRORS ====\n\n" + "\n".join(all_errors))
    if all_warnings:
        with open("mol_creation_warning.log", "w", encoding="utf-8") as fh_warn:
            fh_warn.write("==== MOL CREATION WARNINGS ====\n\n" + "\n".join(all_warnings))

    # ── Summary to console ─────────────────────────────────────────────
    print(f"\n{ANSI_GREEN}Generated {saved_files} MOL files in '{output_dir}'.{ANSI_RESET}")
    print(f"{ANSI_GREEN}Failed to generate {len(all_errors)} MOL files.{ANSI_RESET}")
    if all_errors:
        print(f"{ANSI_RED}See 'mol_creation_error.log' for details.{ANSI_RESET}")
    if all_warnings:
        print(f"{ANSI_ORANGE}See 'mol_creation_warning.log' for fallbacks.{ANSI_RESET}")

    return output_dir
