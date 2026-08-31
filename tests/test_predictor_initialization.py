from __future__ import annotations

import contextlib
import io
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from demiurge_bin import predictor


class PredictorInitializationTests(unittest.TestCase):
    def tearDown(self) -> None:
        predictor.shutdown_persistent_java_processors("test-teardown")

    def test_artifact_validation_exception_is_visible_and_preserved(self):
        failure = RuntimeError("scientific predictor artifact hash mismatch")
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(os.environ, {
            predictor.JAVA_BUILD_DIR_ENV: str(Path(directory) / "build"),
            predictor.JAVA_LIFECYCLE_ENV: "per-batch",
        }, clear=False), mock.patch(
            "demiurge_bin.predictor.verify_predictor_artifacts", side_effect=failure
        ), contextlib.redirect_stderr(io.StringIO()) as stderr:
            with self.assertRaises(RuntimeError) as raised:
                predictor.run_java_batch_processor(
                    Path(directory) / "mols", "1H", 1, "4G", Path(directory) / "raw"
                )
        self.assertIs(raised.exception, failure)
        self.assertIn("artifact hash mismatch", stderr.getvalue())

    def test_javac_failure_is_visible_and_preserved(self):
        failure = subprocess.CalledProcessError(2, ["javac"])
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(os.environ, {
            predictor.JAVA_BUILD_DIR_ENV: str(Path(directory) / "build"),
            predictor.JAVA_LIFECYCLE_ENV: "per-batch",
        }, clear=False), mock.patch(
            "demiurge_bin.predictor.verify_predictor_artifacts", return_value={}
        ), mock.patch(
            "demiurge_bin.predictor._resolve_java_tool", return_value="javac"
        ), mock.patch(
            "demiurge_bin.predictor.subprocess.run", side_effect=failure
        ), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(subprocess.CalledProcessError) as raised:
                predictor.run_java_batch_processor(
                    Path(directory) / "mols", "13C", 1, "4G", Path(directory) / "raw"
                )
        self.assertIs(raised.exception, failure)

    def test_requested_output_directory_is_used_exactly(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(os.environ, {
            predictor.JAVA_LIFECYCLE_ENV: "per-batch",
        }, clear=False), mock.patch(
            "demiurge_bin.predictor._ensure_java_compiled", return_value=("classpath", "predictor.Main")
        ), mock.patch(
            "demiurge_bin.predictor._resolve_java_tool", return_value="java"
        ), mock.patch(
            "demiurge_bin.predictor.subprocess.run",
            return_value=subprocess.CompletedProcess(["java"], 0),
        ) as run:
            requested = Path(directory) / "exact" / "raw_1h"
            result = predictor.run_java_batch_processor(
                Path(directory) / "mols", "1H", 1, "4G", requested
            )
            self.assertEqual(result, str(requested.resolve()))
            self.assertTrue(requested.is_dir())
            command = run.call_args.args[0]
            self.assertEqual(command[command.index(str(Path(directory) / "mols")) + 1], str(requested.resolve()))

    def test_linux_build_cache_is_writable_and_source_hash_deterministic(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(os.environ, {
            predictor.JAVA_BUILD_DIR_ENV: str(Path(directory) / "build"),
        }, clear=False), mock.patch(
            "demiurge_bin.predictor.verify_predictor_artifacts", return_value={}
        ), mock.patch(
            "demiurge_bin.predictor.platform.system", return_value="Linux"
        ), mock.patch(
            "demiurge_bin.predictor._resolve_java_tool", return_value="javac"
        ):
            calls = []

            def compile_once(command, **_kwargs):
                calls.append(command)
                build = Path(command[command.index("-d") + 1])
                (build / "predictor").mkdir(parents=True, exist_ok=True)
                (build / "predictor" / "BatchProcessor1H.class").write_bytes(b"compiled")
                return subprocess.CompletedProcess(command, 0)

            with mock.patch("demiurge_bin.predictor.subprocess.run", side_effect=compile_once):
                first = predictor._ensure_java_compiled("1H")
                second = predictor._ensure_java_compiled("1H")
            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)
            self.assertIn(":", first[0])
            self.assertTrue((Path(directory) / "build" / "predictor" / "BatchProcessor1H.source.sha256").is_file())

    def test_dgx_default_build_directory_is_job_process_isolated(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.dict(os.environ, {
            "SLURM_TMPDIR": directory,
            "SLURM_JOB_ID": "325963",
            "SLURM_ARRAY_TASK_ID": "7",
        }, clear=False):
            os.environ.pop(predictor.JAVA_BUILD_DIR_ENV, None)
            previous = predictor._DEFAULT_BUILD_DIR
            predictor._DEFAULT_BUILD_DIR = None
            try:
                build = predictor._get_build_dir()
                self.assertTrue(build.is_dir())
                self.assertIn("j325963_t7_p", build.name)
                self.assertEqual(build.parent, Path(directory).resolve())
            finally:
                predictor._DEFAULT_BUILD_DIR = previous


if __name__ == "__main__":
    unittest.main()
