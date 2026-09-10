from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from demiurge_bin import bucketing, contracts, preparation
from demiurge_bin.pipeline import _create_ecfp_generator, _ecfp


_SCREEN_ROOT_VALUE = os.environ.get("DEMIURGE_SCREEN_REFERENCE_ROOT", "").strip()
SCREEN_ROOT = Path(_SCREEN_ROOT_VALUE) if _SCREEN_ROOT_VALUE else None


class ScientificContractTests(unittest.TestCase):
    def test_dimensions_and_demiurge_specific_total_order(self):
        self.assertEqual(contracts.MODE_FEATURE_DIMENSIONS, {
            "1H": 200, "13C": 200, "hybrid": 400, "FP": 2048, "total": 2448,
        })
        total = contracts.scientific_contract("total")
        self.assertEqual(total["contract_id"], "DEMIURGE_TOTAL_NMR_V2_H_C_ECFP4")
        self.assertEqual(total["feature_order"], [
            "0:200 1H count buckets",
            "200:400 13C count buckets",
            "400:2448 ECFP4 bits",
        ])

    def test_active_preparation_contains_no_legacy_3d_or_openbabel(self):
        text = Path("demiurge_bin/preparation.py").read_text(encoding="utf-8").lower()
        self.assertNotIn("embedmolecule", text)
        self.assertNotIn("rdcoordgen", text)
        self.assertNotIn("subprocess", text)
        self.assertNotIn("obabel", text)

    def test_predictor_artifacts_are_hash_pinned(self):
        actual = contracts.verify_predictor_artifacts(Path.cwd())
        self.assertEqual(actual, contracts.PREDICTOR_ARTIFACT_SHA256)

    @unittest.skipUnless(
        SCREEN_ROOT is not None and (SCREEN_ROOT / "engine/training_contract.py").is_file(),
        "optional screen reference checkout not configured",
    )
    def test_preparation_is_byte_exact_with_screen_reference(self):
        assert SCREEN_ROOT is not None
        spec = importlib.util.spec_from_file_location(
            "screen_training_contract_reference",
            SCREEN_ROOT / "engine" / "training_contract.py",
        )
        module = importlib.util.module_from_spec(spec)
        assert spec and spec.loader
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        for smiles in ("CCO", "C[C@H](O)Cl", "[NH4+]", "c1ccccc1", "CC(=O)[O-].[Na+]"):
            expected = module.prepare_mol_v2(smiles)
            actual = preparation.prepare_mol_v2(smiles)
            self.assertEqual(actual[0], expected[0])
            self.assertEqual(actual[1].encode("utf-8"), expected[1].encode("utf-8"))
        for invalid in ("not-smiles",):
            with self.assertRaises(Exception) as left:
                module.prepare_mol_v2(invalid)
            with self.assertRaises(type(left.exception)):
                preparation.prepare_mol_v2(invalid)

    def test_bucket_boundaries_are_exact_and_maximum_is_inclusive(self):
        vector, inside, outside = bucketing.bucket_shifts(
            [-1.01, -1.0, 17.0, 17.01], -1.0, 17.0, 200
        )
        self.assertEqual((inside, outside), (2, 2))
        self.assertEqual(float(vector[0]), 1.0)
        self.assertEqual(float(vector[-1]), 1.0)
        self.assertEqual(float(np.sum(vector)), 2.0)

    def test_indexed_prediction_parser_and_nmr_vector(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h = root / "h.csv"
            c = root / "c.csv"
            h.write_text("mol_atom_index,element,parent_atom_index,shift\n1,H,2,1.00\n", encoding="utf-8")
            c.write_text("mol_atom_index,element,shift\n2,C,50.00\n", encoding="utf-8")
            hv, cv, combined, diagnostics = bucketing.vectors_from_raw(h, c)
            self.assertEqual((len(hv), len(cv), len(combined)), (200, 200, 400))
            self.assertEqual(int(np.sum(combined)), 2)
            self.assertEqual(diagnostics["h_in_range"], 1)
            self.assertEqual(diagnostics["c_in_range"], 1)

    def test_ecfp4_is_exactly_the_existing_demiurge_contract(self):
        from rdkit import Chem
        from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

        production = _create_ecfp_generator()
        historical = GetMorganGenerator(radius=2, fpSize=2048)
        for smiles in ("CCO", "c1ccccc1", "N[C@@H](C)C(=O)O", "C[N+](C)(C)C"):
            expected = [int(value) for value in historical.GetFingerprint(Chem.MolFromSmiles(smiles)).ToBitString()]
            self.assertEqual(_ecfp(smiles, production), expected)


if __name__ == "__main__":
    unittest.main()
