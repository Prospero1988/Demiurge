from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from demiurge_bin import preflight


class DeploymentPreflightTests(unittest.TestCase):
    def test_missing_java_is_actionable(self):
        with mock.patch("demiurge_bin.preflight._resolve_java_tool", side_effect=FileNotFoundError("java unavailable")):
            with self.assertRaisesRegex(FileNotFoundError, "java unavailable"):
                preflight.run_preflight()

    def test_missing_javac_is_actionable(self):
        def resolve(name):
            if name == "javac":
                raise FileNotFoundError("javac unavailable")
            return "/conda/bin/java"

        with mock.patch("demiurge_bin.preflight._resolve_java_tool", side_effect=resolve):
            with self.assertRaisesRegex(FileNotFoundError, "javac unavailable"):
                preflight.run_preflight()

    def test_complete_toolchain_preflight_reports_paths_versions_artifacts_and_build(self):
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory) / "build"
            completed = subprocess.CompletedProcess(["tool", "-version"], 0, stdout="", stderr="openjdk 23.0.2\n")
            with mock.patch("demiurge_bin.preflight._resolve_java_tool", side_effect=["/conda/bin/java", "/conda/bin/javac"]), mock.patch(
                "demiurge_bin.preflight.verify_predictor_artifacts", return_value={"predictorh.jar": "abc"}
            ), mock.patch("demiurge_bin.preflight._get_build_dir", return_value=build), mock.patch(
                "demiurge_bin.preflight.subprocess.run", return_value=completed
            ):
                report = preflight.run_preflight(Path(directory))
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["java"]["path"], "/conda/bin/java")
            self.assertEqual(report["javac"]["path"], "/conda/bin/javac")
            self.assertEqual(report["java"]["version"], "openjdk 23.0.2")
            self.assertTrue(report["java_build_directory_writable"])

    def test_conda_environment_pins_full_openjdk(self):
        environment = (Path(__file__).resolve().parents[1] / "conda_environment.yml").read_text(encoding="utf-8")
        self.assertIn("openjdk=23.0.2", environment)
        self.assertNotIn("Java is intentionally external", environment)


if __name__ == "__main__":
    unittest.main()
