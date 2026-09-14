from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from demiurge_bin.pipeline import RunConfig, _create_ecfp_generator, _ecfp, run_pipeline
from demiurge_bin.preparation import NO_SCAFFOLD, molecular_identity_v2, prepare_mol_v2


def fake_java_factory(observed: dict[str, int]):
    def fake_java(mol_directory, nucleus, java_threads=2, java_heap="4G", output_directory=None):
        del java_threads, java_heap
        mols = sorted(Path(mol_directory).glob("*.mol"))
        observed[nucleus] = observed.get(nucleus, 0) + len(mols)
        output = Path(output_directory)
        output.mkdir(parents=True, exist_ok=True)
        for mol in mols:
            text = (
                "mol_atom_index,element,parent_atom_index,shift\n1,H,2,1.00\n"
                if nucleus == "1H" else
                "mol_atom_index,element,shift\n2,C,50.00\n"
            )
            (output / f"{mol.stem}.csv").write_text(text, encoding="utf-8")
        (output / ".spectraprints_unified_profile_branches.jsonl").write_text("", encoding="utf-8")
        return str(output)
    return fake_java


class RecordIdentityTests(unittest.TestCase):
    def test_identity_equivalences_and_distinctions_follow_nmr_v2(self):
        equivalent = [
            ("c1ccccc1", "C1=CC=CC=C1"),
            ("CCO", "OCC"),
            ("CCO", "[CH3][CH2][OH]"),
            ("C[C@H](O)F", "C[C@@H](O)F"),
            ("CCO.[Na+]", "[Na+].OCC"),
        ]
        for left, right in equivalent:
            self.assertEqual(
                molecular_identity_v2(left).identity_sha256,
                molecular_identity_v2(right).identity_sha256,
            )
            self.assertEqual(prepare_mol_v2(left)[1], prepare_mol_v2(right)[1])
            generator = _create_ecfp_generator()
            self.assertEqual(_ecfp(left, generator), _ecfp(right, generator))
        distinct = [
            ("C[N+](C)(C)C", "CN(C)C"),
            ("O=C1NC=CC=C1", "OC1=NC=CC=C1"),
            ("CCO", "CCN"),
            ("CCO", "CCO.[Na+]"),
        ]
        for left, right in distinct:
            self.assertNotEqual(
                molecular_identity_v2(left).identity_sha256,
                molecular_identity_v2(right).identity_sha256,
            )
        self.assertEqual(molecular_identity_v2("CCO").murcko_smiles, NO_SCAFFOLD)

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    def test_compute_dedup_preserves_every_experimental_record(self, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "duplicates.csv"
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(["MOLECULE_NAME", "SMILES", "VALUE"])
                writer.writerows([
                    ("same-name", "c1ccccc1", 5.0),
                    ("same-name", "C1=CC=CC=C1", 5.1),
                    ("same-name", "c1ccncc1", 6.0),
                    ("same-name", "c1ccccc1", 5.0),
                    ("ethanol-a", "CCO", 1.0),
                    ("ethanol-b", "OCC", 2.0),
                ])
            observed: dict[str, int] = {}
            config = RunConfig(
                input_path=source, mode="total", output_root=root / "out",
                temp_root=root / "scratch", label_column=3, prep_workers=1,
                java_threads=2, java_heap="4G", batch_size=2,
                java_lifecycle="persistent", include_murcko=True,
            )
            with mock.patch(
                "demiurge_bin.pipeline.predictor.run_java_batch_processor",
                side_effect=fake_java_factory(observed),
            ):
                summary = run_pipeline(config)
            self.assertEqual((summary["successful"], summary["failed"]), (6, 0))
            self.assertEqual(observed, {"1H": 3, "13C": 3})
            self.assertEqual(summary["record_audit"], {
                "accepted_output_records": 6,
                "unique_molecule_names": 3,
                "duplicate_molecule_name_records": 3,
                "unique_molecular_structures": 3,
                "reused_molecular_structure_records": 3,
                "exact_duplicate_input_records": 1,
                "unique_feature_computations": 3,
            })
            with Path(summary["final_output"]).open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 6)
            self.assertEqual(len({row["RECORD_ID"] for row in rows}), 6)
            self.assertEqual([float(row["LABEL"]) for row in rows], [5.0, 5.1, 6.0, 5.0, 1.0, 2.0])
            self.assertEqual(rows[0]["MURCKO_ID"], rows[1]["MURCKO_ID"])
            self.assertEqual(rows[0]["MURCKO_ID"], rows[3]["MURCKO_ID"])
            self.assertNotEqual(rows[0]["MURCKO_ID"], rows[2]["MURCKO_ID"])
            features = [[row[f"FEATURE_{index}"] for index in range(1, 2449)] for row in rows]
            self.assertEqual(features[0], features[1])
            self.assertEqual(features[0], features[3])
            self.assertEqual(features[4], features[5])
            self.assertNotEqual(features[0], features[2])


if __name__ == "__main__":
    unittest.main()
