from __future__ import annotations

import os
import shutil
import tempfile
import unittest
from pathlib import Path

from demiurge_bin import predictor
from demiurge_bin.bucketing import vectors_from_raw
from demiurge_bin.preparation import prepare_mol_v2


def command_available(name: str) -> bool:
    return shutil.which(name) is not None


@unittest.skipUnless(command_available("java") and command_available("javac"), "JDK unavailable")
class JavaIntegrationTests(unittest.TestCase):
    def test_lifecycle_threads_and_repetitions_are_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mols = root / "mols"
            mols.mkdir()
            for molecule_id, smiles in (("ethanol", "CCO"), ("benzene", "c1ccccc1")):
                _, block = prepare_mol_v2(smiles)
                (mols / f"{molecule_id}.mol").write_text(block, encoding="utf-8", newline="\n")
            os.environ[predictor.JAVA_DIAGNOSTICS_DIR_ENV] = str(root / "diagnostics")
            os.environ[predictor.JAVA_PREDICTOR_MODE_ENV] = "thread-local"
            os.environ["SPECTRAPRINTS_UNIFIED_PROFILE"] = "1"
            outputs = []
            try:
                os.environ[predictor.JAVA_LIFECYCLE_ENV] = "per-batch"
                baseline_h = root / "baseline_h"; baseline_c = root / "baseline_c"
                self.assertIsNotNone(predictor.run_java_batch_processor(mols, "1H", 1, "4G", baseline_h))
                self.assertIsNotNone(predictor.run_java_batch_processor(mols, "13C", 1, "4G", baseline_c))
                outputs.append((baseline_h, baseline_c))
                os.environ[predictor.JAVA_LIFECYCLE_ENV] = "persistent"
                for label, threads in (("2t_first", 2), ("2t_repeat", 2), ("4t", 4)):
                    h = root / f"persistent_{label}_h"
                    c = root / f"persistent_{label}_c"
                    self.assertIsNotNone(predictor.run_java_batch_processor(mols, "1H", threads, "4G", h))
                    self.assertIsNotNone(predictor.run_java_batch_processor(mols, "13C", threads, "4G", c))
                    outputs.append((h, c))
            finally:
                predictor.shutdown_persistent_java_processors("test-finally")
            reference_h, reference_c = outputs[0]
            for h_dir, c_dir in outputs[1:]:
                for molecule_id in ("ethanol", "benzene"):
                    self.assertEqual(
                        (reference_h / f"{molecule_id}.csv").read_bytes(),
                        (h_dir / f"{molecule_id}.csv").read_bytes(),
                    )
                    self.assertEqual(
                        (reference_c / f"{molecule_id}.csv").read_bytes(),
                        (c_dir / f"{molecule_id}.csv").read_bytes(),
                    )
                    reference = vectors_from_raw(
                        reference_h / f"{molecule_id}.csv",
                        reference_c / f"{molecule_id}.csv",
                    )[2]
                    candidate = vectors_from_raw(
                        h_dir / f"{molecule_id}.csv",
                        c_dir / f"{molecule_id}.csv",
                    )[2]
                    self.assertEqual(reference.tobytes(), candidate.tobytes())


if __name__ == "__main__":
    unittest.main()
