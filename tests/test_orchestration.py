from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from demiurge_bin.run_state import atomic_write_json
from demiurge_supervisor import (
    MAX_ATTEMPTS,
    campaign_status,
    classify_failure,
    ensure_previous_arrays_inactive,
    load_manifest,
    prepare_manifest,
    retry_decision,
    submit_chain,
)
from orchestration.staging import cleanup_staging, stage_input


class OrchestrationTests(unittest.TestCase):
    def args(self, root: Path):
        return SimpleNamespace(
            project_root=Path.cwd(), input_dir=root / "inputs", pattern="*.csv",
            output_root=root / "campaigns", scratch_root=root / "scratch",
            campaign="test_campaign", mode="total", label_column=3,
            batch_size=1000, prep_workers=4, java_threads=2, java_heap="4G",
            java_lifecycle="persistent", cpus_per_task=6, memory="16G",
            partition="dgx_long", time="24:00:00", max_concurrent=3,
            job_name_prefix="DEMIURGE test", conda_root="/raid/soft/miniconda",
            conda_env="demiurge", no_staging=False,
            retain_scientific_artifacts=False, dry_run=True,
        )

    def create_inputs(self, root: Path):
        directory = root / "inputs"
        directory.mkdir()
        for name in ("b.csv", "a.csv"):
            (directory / name).write_text("MOLECULE_NAME,SMILES,LABEL\nx,CCO,1\n", encoding="utf-8")

    def test_manifest_discovery_is_sorted_and_three_attempt_chain_is_finite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.create_inputs(root)
            path, manifest = prepare_manifest(self.args(root))
            self.assertEqual([Path(task["input_path"]).name for task in manifest["tasks"]], ["a.csv", "b.csv"])
            jobs = submit_chain(path, manifest, [0, 1], dry_run=True)
            self.assertEqual(jobs, ["DRYRUN1", "DRYRUN2", "DRYRUN3"])
            self.assertEqual(len(manifest["submissions"]), MAX_ATTEMPTS)
            self.assertTrue((path.parent / "logs").is_dir())
            self.assertEqual(manifest["submissions"][1]["dependency"], "DRYRUN1")
            _, loaded = load_manifest(path)
            self.assertEqual(loaded["max_attempts"], 3)

    def test_retry_classification_is_fail_closed(self):
        self.assertEqual(classify_failure("broken pipe to Java"), ("TRANSIENT_OPERATIONAL", True))
        self.assertEqual(classify_failure("scientific configuration differs"), ("PERMANENT_CONFIGURATION", False))
        self.assertEqual(classify_failure("mysterious failure"), ("UNCLASSIFIED_FAIL_CLOSED", False))

    def test_retry_decision_stops_done_and_permanent_tasks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.create_inputs(root)
            path, manifest = prepare_manifest(self.args(root))
            output = Path(manifest["tasks"][0]["output_root"])
            output.mkdir(parents=True)
            atomic_write_json(output / "checkpoint.json", {"status": "DONE"})
            self.assertEqual(retry_decision(path, manifest, 0, 2), "SKIP_DONE")
            output2 = Path(manifest["tasks"][1]["output_root"])
            output2.mkdir(parents=True)
            atomic_write_json(output2 / "checkpoint.json", {"status": "FAILED"})
            attempt = path.parent / "attempt_history" / "task-0001" / "attempt-01.json"
            atomic_write_json(attempt, {"retryable": False})
            self.assertEqual(retry_decision(path, manifest, 1, 2), "SKIP_PERMANENT")

    def test_historical_invalid_job_is_inactive_but_active_and_unknown_errors_block(self):
        manifest = {"submissions": [{"job_id": "123", "attempt": 1}]}

        def missing(command, **_kwargs):
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="slurm_load_jobs error: Invalid job id specified")

        ensure_previous_arrays_inactive(manifest, missing)

        def active(command, **_kwargs):
            return subprocess.CompletedProcess(command, 0, stdout="123|RUNNING\n", stderr="")

        with self.assertRaisesRegex(RuntimeError, "still active"):
            ensure_previous_arrays_inactive(manifest, active)

        def unexpected(command, **_kwargs):
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="controller unavailable")

        with self.assertRaisesRegex(RuntimeError, "Could not verify"):
            ensure_previous_arrays_inactive(manifest, unexpected)

    def test_campaign_status_uses_manifest_and_checkpoints(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.create_inputs(root)
            path, manifest = prepare_manifest(self.args(root))
            output = Path(manifest["tasks"][0]["output_root"])
            output.mkdir(parents=True)
            atomic_write_json(output / "checkpoint.json", {"status": "DONE", "successful": 1, "failed": 0})
            text, progress = campaign_status(path, manifest)
            self.assertIn("total_expected=2", text)
            self.assertIn("completed_shards=1", text)
            self.assertIn("pending_shards=1", text)
            self.assertTrue(progress.is_file())

    def test_staging_copy_hash_and_marker_protected_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            source.write_text("a,b\n1,2\n", encoding="utf-8")
            result = stage_input(source, root / "scratch", "campaign", "1", 0, 1)
            staged = Path(result["runtime_input"])
            self.assertEqual(staged.read_bytes(), source.read_bytes())
            with self.assertRaisesRegex(RuntimeError, "marker mismatch"):
                cleanup_staging(root / "scratch", Path(result["staging_directory"]), "wrong")
            cleanup_staging(root / "scratch", Path(result["staging_directory"]), result["owner_token"])
            self.assertFalse(Path(result["staging_directory"]).exists())

    def test_shell_worker_is_real_bash_and_contains_no_wrap(self):
        text = Path("orchestration/slurm_worker.sh").read_text(encoding="utf-8")
        self.assertTrue(text.startswith("#!/usr/bin/env bash\n"))
        self.assertIn("set -Eeuo pipefail", text)
        self.assertNotIn("--wrap", text)
        self.assertIn("python demiurge.py", text)


if __name__ == "__main__":
    unittest.main()
