from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from demiurge_bin.pipeline import RunConfig, read_status, resume_pipeline, run_pipeline


def write_input(path: Path, count: int = 2, invalid: bool = False) -> None:
    rows = [(f"mol-{index}", "not-smiles" if invalid and index == count - 1 else "CCO", 5.0 + index) for index in range(count)]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["MOLECULE_NAME", "SMILES", "VALUE"])
        writer.writerows(rows)


def fake_java(mol_directory, nucleus, java_threads=8, java_heap="4G", output_directory=None):
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    for mol in sorted(Path(mol_directory).glob("*.mol")):
        if nucleus == "1H":
            text = "mol_atom_index,element,parent_atom_index,shift\n1,H,2,1.00\n"
        else:
            text = "mol_atom_index,element,shift\n2,C,50.00\n"
        (output / f"{mol.stem}.csv").write_text(text, encoding="utf-8")
    (output / ".spectraprints_unified_profile_branches.jsonl").write_text("", encoding="utf-8")
    return str(output)


class PipelineTests(unittest.TestCase):
    def config(self, root: Path, input_path: Path, backend: str = "local", **values):
        options = dict(
            input_path=input_path,
            mode="total",
            output_root=root,
            temp_root=root.parent / "scratch",
            label_column=3,
            prep_workers=1,
            java_threads=2,
            java_heap="4G",
            batch_size=1,
            java_lifecycle="persistent",
            retain_scientific_artifacts=True,
            backend=backend,
        )
        options.update(values)
        return RunConfig(**options)

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    @mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java)
    def test_total_composition_and_failure_reporting(self, _java, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            write_input(source, 2, invalid=True)
            summary = run_pipeline(self.config(root / "output", source))
            self.assertEqual((summary["successful"], summary["failed"]), (1, 1))
            final = Path(summary["final_output"])
            with final.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(len(rows[0]), 2450)
            self.assertEqual(len(rows), 2)
            self.assertTrue((root / "output" / "failures.jsonl").read_text(encoding="utf-8"))

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    @mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java)
    def test_local_and_slurm_worker_backends_are_byte_exact(self, _java, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            write_input(source, 2)
            local = run_pipeline(self.config(root / "local", source, "local"))
            slurm = run_pipeline(self.config(root / "slurm", source, "slurm-worker"))
            self.assertEqual(list((root / "scratch").iterdir()), [])
            self.assertEqual(Path(local["final_output"]).read_bytes(), Path(slurm["final_output"]).read_bytes())
            for index in range(2):
                left = root / "local" / "batches" / f"batch_{index:08d}"
                right = root / "slurm" / "batches" / f"batch_{index:08d}"
                self.assertEqual((left / "features.csv").read_bytes(), (right / "features.csv").read_bytes())
                self.assertEqual((left / "metadata.jsonl").read_bytes(), (right / "metadata.jsonl").read_bytes())
                for name in ("mols", "raw_1h", "raw_13c"):
                    left_files = sorted((left / "scientific_artifacts" / name).glob("*"))
                    right_files = sorted((right / "scientific_artifacts" / name).glob("*"))
                    self.assertEqual([p.name for p in left_files], [p.name for p in right_files])
                    for lpath, rpath in zip(left_files, right_files):
                        self.assertEqual(lpath.read_bytes(), rpath.read_bytes())
                for name in ("raw_1h", "raw_13c"):
                    branch = ".spectraprints_unified_profile_branches.jsonl"
                    self.assertEqual(
                        (left / "scientific_artifacts" / name / branch).read_bytes(),
                        (right / "scientific_artifacts" / name / branch).read_bytes(),
                    )

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    def test_transient_batch_retry_is_bounded_and_recovers(self, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            write_input(source, 1)
            calls = {"count": 0}

            def flaky(*args, **kwargs):
                calls["count"] += 1
                if calls["count"] <= 2:
                    return None
                return fake_java(*args, **kwargs)

            with mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=flaky):
                summary = run_pipeline(self.config(root / "output", source, mode="1H"))
            self.assertEqual(summary["successful"], 1)
            self.assertEqual(calls["count"], 3)

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    def test_checkpoint_resume_replays_only_uncommitted_batch(self, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            write_input(source, 2)

            def fail_second(mol_directory, nucleus, **kwargs):
                if any(path.stem == "m00000001" for path in Path(mol_directory).glob("*.mol")):
                    return None
                return fake_java(mol_directory, nucleus, **kwargs)

            with mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fail_second):
                with self.assertRaises(RuntimeError):
                    run_pipeline(self.config(root / "output", source, mode="1H"))
            checkpoint = json.loads((root / "output" / "checkpoint.json").read_text(encoding="utf-8"))
            self.assertEqual(checkpoint["next_row_index"], 1)
            self.assertEqual(checkpoint["committed_batches"], 1)
            with mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java):
                summary = resume_pipeline(
                    root / "output", root / "resume-scratch",
                    prep_workers=1, java_threads=2, batch_size=1, java_lifecycle="persistent",
                )
            with Path(summary["final_output"]).open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual([row[0] for row in rows[1:]], ["mol-0", "mol-1"])

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    @mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java)
    def test_status_reads_durable_checkpoint(self, _java, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            write_input(source, 1)
            run_pipeline(self.config(root / "output", source, mode="1H"))
            text, path = read_status(root / "output")
            self.assertIn("status=DONE", text)
            self.assertIn("percentage=100.000", text)
            self.assertEqual(path.name, "production_progress.txt")

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    @mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java)
    def test_resume_rejects_changed_scientific_identity(self, _java, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.csv"
            write_input(source, 1)
            run_pipeline(self.config(root / "output", source, mode="1H"))
            checkpoint = json.loads((root / "output" / "checkpoint.json").read_text(encoding="utf-8"))
            checkpoint["status"] = "FAILED"
            checkpoint["scientific_config"]["contract"]["feature_dimension"] = 201
            (root / "output" / "checkpoint.json").write_text(
                json.dumps(checkpoint), encoding="utf-8"
            )
            with self.assertRaisesRegex(RuntimeError, "Scientific configuration differs"):
                resume_pipeline(root / "output", root / "resume-scratch")


if __name__ == "__main__":
    unittest.main()
