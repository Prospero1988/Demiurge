from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import demiurge_performance_gate as gate


def write_summary(root: Path, *, throughput: float, failed: int = 0) -> None:
    root.mkdir(parents=True)
    (root / "summary.json").write_text(json.dumps({
        "status": "DONE",
        "contract_id": "DEMIURGE_TOTAL_NMR_V2_H_C_ECFP4",
        "total": 10,
        "successful": 10 - failed,
        "failed": failed,
        "wall_time_seconds": 10.0 / throughput,
        "molecules_per_second": throughput,
        "stage_timing_seconds": {"preparation": 1.0, "java_1h": 2.0, "java_13c": 3.0, "features": 0.1},
        "operational": {"java_threads": 2},
    }), encoding="utf-8")


class PerformanceGateTests(unittest.TestCase):
    def test_only_exact_qc_clean_runs_are_ranked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fast, slow, failed, output = root / "fast", root / "slow", root / "failed", root / "aggregate"
            write_summary(fast, throughput=5.0)
            write_summary(slow, throughput=2.0)
            write_summary(failed, throughput=9.0, failed=1)
            argv = ["demiurge_performance_gate.py", "--run", f"fast={fast}", "--run", f"slow={slow}", "--run", f"failed={failed}", "--output-root", str(output)]
            with mock.patch.object(sys, "argv", argv):
                self.assertEqual(gate.main(), 0)
            document = json.loads((output / "performance_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(document["ranking"], ["fast", "slow"])
            with (output / "performance_summary.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertTrue((output / "performance_summary.md").is_file())


if __name__ == "__main__":
    unittest.main()
