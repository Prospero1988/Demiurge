from __future__ import annotations

import csv
import json
import sqlite3
import tempfile
import unittest
from contextlib import closing
from pathlib import Path
from unittest import mock

import numpy as np

from demiurge_bin.io import create_input_reader
from demiurge_bin.io.factory import create_output_writer
from demiurge_bin.contracts import scientific_contract
from demiurge_bin.pipeline import RunConfig, resume_pipeline, run_pipeline


def fake_java(mol_directory, nucleus, java_threads=8, java_heap="4G", output_directory=None):
    del java_threads, java_heap
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


def write_csv(path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["compound_id", "canonical_smiles", "activity"])
        writer.writerows([
            ("mol-a", "CCN", 5.25),
            ("mol-b", "CCO", 6.5),
            ("mol-bad", "not-smiles", 1.0),
        ])


def write_sqlite(path: Path) -> None:
    with closing(sqlite3.connect(path)) as connection:
        connection.execute("CREATE TABLE compounds(compound_id TEXT, canonical_smiles TEXT, activity)")
        connection.executemany(
            "INSERT INTO compounds VALUES (?,?,?)",
            [("mol-a", "CCN", 5.25), ("mol-b", "CCO", 6.5), ("mol-bad", "not-smiles", 1.0)],
        )
        connection.commit()


def csv_rows(path: Path) -> dict[str, tuple[object, bytes]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader)
        return {
            row[0]: (float(row[1]), np.asarray(row[2:], dtype="<f4").tobytes())
            for row in reader
        }


def sqlite_rows(path: Path, table: str = "features") -> dict[str, tuple[object, bytes]]:
    with closing(sqlite3.connect(path)) as connection:
        return {
            molecule_id: (label, bytes(blob))
            for molecule_id, label, blob in connection.execute(
                f'SELECT molecule_id,label,feature_blob FROM "{table}" WHERE status="SUCCESS"'
            )
        }


class GenericIoTests(unittest.TestCase):
    def config(self, root: Path, source: Path, *, input_format: str, output_format: str) -> RunConfig:
        return RunConfig(
            input_path=source,
            input_format=input_format,
            input_table="compounds" if input_format == "sqlite" else None,
            id_column="compound_id",
            smiles_column="canonical_smiles",
            label_column="activity",
            mode="total",
            output_root=root,
            output_format=output_format,
            output_db=(root / "result.db") if output_format == "sqlite" else None,
            output_table="features" if output_format == "sqlite" else "demiurge_features",
            metadata_table="schema_info" if output_format == "sqlite" else "demiurge_metadata",
            temp_root=root.parent / "scratch",
            prep_workers=1,
            java_threads=2,
            java_heap="4G",
            batch_size=2,
            java_lifecycle="persistent",
        )

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    @mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java)
    def test_all_four_format_combinations_are_scientifically_exact(self, _java, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_csv = root / "molecules.csv"
            source_db = root / "molecules.db"
            write_csv(source_csv)
            write_sqlite(source_db)

            csv_csv = run_pipeline(self.config(root / "csv_csv", source_csv, input_format="csv", output_format="csv"))
            csv_db = run_pipeline(self.config(root / "csv_db", source_csv, input_format="csv", output_format="sqlite"))
            db_csv = run_pipeline(self.config(root / "db_csv", source_db, input_format="sqlite", output_format="csv"))
            db_db = run_pipeline(self.config(root / "db_db", source_db, input_format="sqlite", output_format="sqlite"))

            expected = csv_rows(Path(csv_csv["final_output"]))
            self.assertEqual(expected, sqlite_rows(Path(csv_db["final_output"])))
            self.assertEqual(expected, csv_rows(Path(db_csv["final_output"])))
            self.assertEqual(expected, sqlite_rows(Path(db_db["final_output"])))
            self.assertEqual(Path(csv_csv["final_output"]).read_bytes(), Path(db_csv["final_output"]).read_bytes())
            self.assertEqual(set(expected), {"mol-a", "mol-b"})

            with closing(sqlite3.connect(csv_db["final_output"])) as connection:
                failed = connection.execute(
                    'SELECT molecule_id,label,status,error,feature_blob FROM "features" WHERE status="FAILED"'
                ).fetchone()
                metadata = dict(connection.execute('SELECT key,value_json FROM "schema_info"'))
            self.assertEqual(failed[:3], ("mol-bad", 1.0, "FAILED"))
            self.assertIn("RDKit rejected SMILES", failed[3])
            self.assertIsNone(failed[4])
            self.assertEqual(json.loads(metadata["feature_count"]), 2448)
            self.assertEqual(json.loads(metadata["feature_blob_bytes"]), 9792)
            self.assertEqual(json.loads(metadata["dtype"]), "<f4")
            self.assertEqual(json.loads(metadata["feature_order"]), [
                "0:200 1H count buckets", "200:400 13C count buckets", "400:2448 ECFP4 bits",
            ])
            self.assertEqual(json.loads(metadata["feature_offsets"]), {
                "1H": [0, 200], "13C": [200, 400], "ECFP4": [400, 2448],
            })
            for summary in (csv_csv, csv_db, db_csv, db_db):
                self.assertEqual((summary["successful"], summary["failed"]), (2, 1))

    def test_sqlite_query_and_bounded_batches(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "molecules.db"
            write_sqlite(source)
            reader = create_input_reader(
                source,
                input_format="sqlite",
                input_table=None,
                input_query="SELECT compound_id, canonical_smiles, activity FROM compounds WHERE activity > 1",
                id_column="compound_id",
                smiles_column="canonical_smiles",
                label_column="activity",
            )
            self.assertEqual(reader.describe().total_rows, 2)
            batches = list(reader.iter_batches(1))
            self.assertEqual([len(records) for _, records in batches], [1, 1])
            self.assertEqual([records[0]["molecule_name"] for _, records in batches], ["mol-a", "mol-b"])

    def test_sqlite_configuration_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "molecules.db"
            write_sqlite(source)
            with self.assertRaisesRegex(ValueError, "exactly one"):
                create_input_reader(
                    source,
                    input_format="sqlite",
                    input_table="compounds",
                    input_query="SELECT * FROM compounds",
                    id_column="compound_id",
                    smiles_column="canonical_smiles",
                    label_column="activity",
                )
            (root / "molecules.db-wal").write_bytes(b"active")
            reader = create_input_reader(
                source,
                input_format="sqlite",
                input_table="compounds",
                input_query=None,
                id_column="compound_id",
                smiles_column="canonical_smiles",
                label_column="activity",
            )
            with self.assertRaisesRegex(RuntimeError, "active -wal sidecar"):
                reader.describe()

    def test_existing_sqlite_output_is_not_overwritten_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "existing.db"
            original = b"do-not-replace"
            output.write_bytes(original)
            writer = create_output_writer(
                output_format="sqlite",
                output_root=root / "run",
                output_db=output,
                output_table="features",
                metadata_table="metadata",
                input_stem="input",
                mode="total",
                feature_contract=scientific_contract("total"),
                overwrite_output=False,
            )
            with self.assertRaises(FileExistsError):
                writer.initialize(resume=False)
            self.assertEqual(output.read_bytes(), original)

    def test_input_database_cannot_be_selected_as_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "molecules.db"
            write_sqlite(source)
            config = self.config(Path(directory) / "run", source, input_format="sqlite", output_format="sqlite")
            config = RunConfig(**{**config.__dict__, "output_db": source, "overwrite_output": True})
            with self.assertRaisesRegex(ValueError, "must not be the input"):
                config.validated()

    @mock.patch("demiurge_bin.pipeline.predictor.shutdown_persistent_java_processors")
    @mock.patch("demiurge_bin.pipeline.predictor.run_java_batch_processor", side_effect=fake_java)
    def test_sqlite_resume_replays_only_uncommitted_rows(self, java, _shutdown):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "molecules.csv"
            with source.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(["compound_id", "canonical_smiles", "activity"])
                writer.writerows([
                    ("mol-a", "CCO", 1.0),
                    ("mol-b", "CCN", 2.0),
                    ("mol-c", "CCC", 3.0),
                ])
            config = self.config(root / "output", source, input_format="csv", output_format="sqlite")

            original = java.side_effect
            calls = {"count": 0}

            def fail_second(*args, **kwargs):
                calls["count"] += 1
                if any(path.stem == "m00000002" for path in Path(args[0]).glob("*.mol")):
                    raise RuntimeError("transient test failure")
                return original(*args, **kwargs)

            java.side_effect = fail_second
            with self.assertRaises(RuntimeError):
                run_pipeline(config)
            java.side_effect = fake_java
            summary = resume_pipeline(
                root / "output",
                root / "resume-scratch",
                prep_workers=1,
                java_threads=2,
                batch_size=2,
                java_lifecycle="persistent",
            )
            self.assertEqual((summary["successful"], summary["failed"]), (3, 0))
            self.assertEqual(set(sqlite_rows(Path(summary["final_output"]))), {"mol-a", "mol-b", "mol-c"})


if __name__ == "__main__":
    unittest.main()
