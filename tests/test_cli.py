from __future__ import annotations

import unittest

import demiurge
import demiurge_supervisor


class CliTests(unittest.TestCase):
    def test_existing_csv_run_defaults_remain_compatible(self):
        args = demiurge.build_parser().parse_args([
            "run", "--input", "input.csv", "--mode", "total", "--output-root", "results",
        ])
        self.assertEqual(args.input_format, "csv")
        self.assertEqual(args.output_format, "csv")
        self.assertEqual(args.id_column, "MOLECULE_NAME")
        self.assertEqual(args.smiles_column, "SMILES")
        self.assertEqual(args.label_column, 3)
        self.assertEqual((args.java_threads, args.prep_workers, args.java_heap), (2, 4, "4G"))
        self.assertEqual(args.java_lifecycle, "persistent")

    def test_sqlite_column_names_and_output_options_parse(self):
        args = demiurge.build_parser().parse_args([
            "run", "--input", "input.db", "--input-format", "sqlite",
            "--input-table", "compounds", "--id-column", "compound_id",
            "--smiles-column", "canonical_smiles", "--label-column", "activity",
            "--mode", "total", "--output-root", "results", "--output-format", "sqlite",
            "--output-db", "results/features.db", "--output-table", "vectors",
            "--metadata-table", "vector_schema",
        ])
        self.assertEqual(args.label_column, "activity")
        self.assertEqual(args.output_table, "vectors")
        self.assertEqual(args.metadata_table, "vector_schema")

    def test_supervisor_keeps_csv_defaults_and_accepts_sqlite(self):
        parser = demiurge_supervisor.build_parser()
        base = [
            "submit", "--project-root", ".", "--input-dir", "inputs",
            "--output-root", "outputs", "--scratch-root", "scratch", "--campaign", "test",
        ]
        csv_args = parser.parse_args(base)
        self.assertEqual((csv_args.pattern, csv_args.input_format, csv_args.output_format), ("*.csv", "csv", "csv"))
        db_args = parser.parse_args(base + [
            "--pattern", "*.db", "--input-format", "sqlite", "--input-table", "molecules",
            "--id-column", "id", "--smiles-column", "smi", "--label-column", "value",
            "--output-format", "sqlite",
        ])
        self.assertEqual((db_args.input_table, db_args.label_column, db_args.output_format), ("molecules", "value", "sqlite"))


if __name__ == "__main__":
    unittest.main()
