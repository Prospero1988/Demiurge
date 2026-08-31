from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

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


if __name__ == "__main__":
    unittest.main()
