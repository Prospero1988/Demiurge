from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import demiurge_nmr_v2_gate as gate


def create_output(root: Path, changed: bool = False) -> None:
    for directory in ("mols", "raw_1h", "raw_13c"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    (root / "mols" / "m1.mol").write_bytes(b"mol\n")
    (root / "raw_1h" / "m1.csv").write_bytes(b"h\n" if not changed else b"changed\n")
    (root / "raw_13c" / "m1.csv").write_bytes(b"c\n")
    for directory in ("raw_1h", "raw_13c"):
        (root / directory / gate.BRANCH_FILE).write_bytes(b"{}\n")
    for name, value in (
        ("canonical_identity.json", {"m1": "CCO"}),
        ("preparation_failures.json", {"failures": []}),
        ("nmr_vectors.json", {"m1": {"H_C": [1, 0]}}),
        ("prediction_failures.json", {"failures": []}),
    ):
        (root / name).write_text(json.dumps(value), encoding="utf-8")


class ParityComparatorTests(unittest.TestCase):
    def test_screen_emission_requires_explicit_reference_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(
                implementation="screen",
                screen_root=None,
                corpus=Path(directory) / "missing.jsonl",
                output_root=Path(directory) / "output",
            )
            with self.assertRaisesRegex(ValueError, "--screen-root is required"):
                gate.emit(args)
            self.assertFalse(args.output_root.exists())

    def test_exact_outputs_pass_and_any_raw_change_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left, right = root / "screen", root / "demiurge"
            create_output(left); create_output(right)
            self.assertEqual(gate.compare(SimpleNamespace(screen_output=left, demiurge_output=right)), 0)
            (right / "raw_1h" / "m1.csv").write_bytes(b"changed\n")
            with self.assertRaisesRegex(RuntimeError, "raw 1H"):
                gate.compare(SimpleNamespace(screen_output=left, demiurge_output=right))

    def test_scientific_json_change_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left, right = root / "screen", root / "demiurge"
            create_output(left); create_output(right)
            (right / "canonical_identity.json").write_text('{"m1":"CCC"}', encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "scientific JSON"):
                gate.compare(SimpleNamespace(screen_output=left, demiurge_output=right))

    def test_frozen_reference_hashes_are_fail_closed(self):
        reference = Path(__file__).resolve().parents[1] / "validation" / "frozen_nmr_v2_expected"
        corpus = Path(__file__).resolve().parents[1] / "validation" / "nmr_v2_parity_corpus.jsonl"
        self.assertEqual(gate.verify_frozen_reference(reference, corpus)["contract"], "SPECTRAPRINTS_NMR_V2")
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / "reference"
            import shutil
            shutil.copytree(reference, copied)
            target = copied / "nmr" / "raw_1h" / "ordinary_ethanol.csv"
            target.write_bytes(target.read_bytes() + b"tampered")
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                gate.verify_frozen_reference(copied, corpus)

    def test_standalone_validation_has_no_screen_checkout_dependency(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference"
            candidate = root / "candidate"
            corpus = root / "corpus.jsonl"
            corpus.write_text('{"molecule_id":"m1","smiles":"CCO"}\n', encoding="utf-8")
            create_output(reference / "nmr")
            inventory = {}
            for path in sorted((reference / "nmr").rglob("*")):
                if path.is_file():
                    inventory[path.relative_to(reference).as_posix()] = gate.file_sha256(path)
            (reference / gate.FROZEN_MANIFEST).write_text(json.dumps({
                "schema_version": 1,
                "contract": "SPECTRAPRINTS_NMR_V2",
                "corpus_sha256": gate.file_sha256(corpus),
                "files_sha256": inventory,
            }), encoding="utf-8")

            def fake_emit(args):
                create_output(args.output_root)
                return 0

            args = SimpleNamespace(
                reference_root=reference, corpus=corpus, output_root=candidate,
                java_threads=2, java_heap="4G", java_lifecycle="persistent",
            )
            with mock.patch.object(gate, "emit", side_effect=fake_emit), mock.patch.object(
                gate, "_screen_modules", side_effect=AssertionError("screen checkout was accessed")
            ):
                self.assertEqual(gate.validate_standalone(args), 0)

    def test_validation_sbatch_requires_no_screen_project_root(self):
        script = (Path(__file__).resolve().parents[1] / "benchmarks" / "demiurge_nmr_v2_validation.sbatch").read_text(encoding="utf-8")
        self.assertNotIn("SCREEN_PROJECT_ROOT", script)
        self.assertIn("validate-standalone", script)
        self.assertIn("frozen_nmr_v2_expected", script)


if __name__ == "__main__":
    unittest.main()
