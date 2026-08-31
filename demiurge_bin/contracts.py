"""Frozen scientific contracts for the production Demiurge representations.

The NMR portion intentionally mirrors the validated SPECTRAPRINTS_NMR_V2
training/screening contract.  Operational settings (backend, paths, batching,
thread counts, heap size and lifecycle) are deliberately absent.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


NMR_REPRESENTATION_VERSION = "SPECTRAPRINTS_NMR_V2"
DEMIURGE_1H_CONTRACT = "DEMIURGE_1H_NMR_V2"
DEMIURGE_13C_CONTRACT = "DEMIURGE_13C_NMR_V2"
DEMIURGE_HYBRID_CONTRACT = "DEMIURGE_HYBRID_NMR_V2_H_C"
DEMIURGE_FP_CONTRACT = "DEMIURGE_ECFP4"
DEMIURGE_TOTAL_CONTRACT = "DEMIURGE_TOTAL_NMR_V2_H_C_ECFP4"

H_MIN = -1.0
H_MAX = 17.0
H_BINS = 200
C_MIN = -10.0
C_MAX = 230.0
C_BINS = 200
ECFP_RADIUS = 2
ECFP_BITS = 2048
ECFP_USE_CHIRALITY = False

MODE_CONTRACTS = {
    "1H": DEMIURGE_1H_CONTRACT,
    "13C": DEMIURGE_13C_CONTRACT,
    "hybrid": DEMIURGE_HYBRID_CONTRACT,
    "FP": DEMIURGE_FP_CONTRACT,
    "total": DEMIURGE_TOTAL_CONTRACT,
}

MODE_FEATURE_DIMENSIONS = {
    "1H": H_BINS,
    "13C": C_BINS,
    "hybrid": H_BINS + C_BINS,
    "FP": ECFP_BITS,
    "total": H_BINS + C_BINS + ECFP_BITS,
}

PREDICTOR_ARTIFACT_SHA256 = {
    "cdk-2.9.jar": "60710218b8f9fd206e6151122e630c281462e9588e4b7a279c49c1532a8aeffe",
    "cdk-builder3d-2.9.jar": "2c3add480bc7363b5fe6da076f873543b47355630149927b9420126af78542ca",
    "predictorc.jar": "e3c3365fb3ffdccd79bb1c39c457c2486e6170f88eeaca5f36c09587950a5090",
    "predictorh.jar": "529e2c89279aaafcf63347460775693d0ff5120d17dae051e31b9e55f6d1e67d",
}


def scientific_contract(mode: str) -> dict[str, Any]:
    if mode not in MODE_CONTRACTS:
        raise ValueError(f"Unsupported Demiurge mode: {mode}")
    feature_order: list[str]
    if mode == "1H":
        feature_order = ["0:200 1H count buckets"]
    elif mode == "13C":
        feature_order = ["0:200 13C count buckets"]
    elif mode == "hybrid":
        feature_order = ["0:200 1H count buckets", "200:400 13C count buckets"]
    elif mode == "FP":
        feature_order = ["0:2048 ECFP4 bits"]
    else:
        feature_order = [
            "0:200 1H count buckets",
            "200:400 13C count buckets",
            "400:2448 ECFP4 bits",
        ]
    return {
        "contract_id": MODE_CONTRACTS[mode],
        "nmr_representation_version": (
            NMR_REPRESENTATION_VERSION if mode != "FP" else None
        ),
        "preparation": {
            "input": "raw SMILES",
            "canonicalization": "RDKit canonical isomeric SMILES",
            "explicit_hydrogens": True,
            "python_etkdg": False,
            "python_3d": False,
            "layout": "RDKit Compute2DCoords",
            "stereochemistry_before_write": "removed",
            "writer": "RDKit V3000",
            "openbabel_fallback": False,
        },
        "java_nmr": {
            "prediction_mode": "3D-first",
            "reconstruction": "CDK ModelBuilder3D",
            "solvent": "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)",
            "predictor_lifecycle": "thread-local PredictionTool; persistent JVM supported",
        },
        "nmr": {
            "aggregation": "unnormalized count buckets",
            "1H": {"minimum_ppm": H_MIN, "maximum_ppm": H_MAX, "bins": H_BINS, "maximum_inclusive": True},
            "13C": {"minimum_ppm": C_MIN, "maximum_ppm": C_MAX, "bins": C_BINS, "maximum_inclusive": True},
        },
        "ecfp": {
            "name": "ECFP4",
            "radius": ECFP_RADIUS,
            "n_bits": ECFP_BITS,
            "use_chirality": ECFP_USE_CHIRALITY,
        },
        "feature_order": feature_order,
        "feature_dimension": MODE_FEATURE_DIMENSIONS[mode],
        "feature_dtype": "float32",
    }


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def object_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_predictor_artifacts(project_root: Path) -> dict[str, str]:
    predictor = project_root / "predictor"
    actual: dict[str, str] = {}
    for filename, expected in PREDICTOR_ARTIFACT_SHA256.items():
        path = predictor / filename
        if not path.is_file():
            raise FileNotFoundError(f"Required scientific predictor artifact is missing: {path}")
        digest = file_sha256(path)
        if digest != expected:
            raise RuntimeError(
                f"Scientific predictor artifact hash mismatch for {path}: "
                f"expected={expected} actual={digest}"
            )
        actual[filename] = digest
    return actual
