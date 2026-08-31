#!/usr/bin/env python3
"""Aggregate completed Demiurge validation runs without rerunning science."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from demiurge_bin.run_state import atomic_write_json, atomic_write_text, utc_now


def load_run(value: str) -> dict:
    if "=" not in value:
        raise ValueError("Each --run must be LABEL=OUTPUT_ROOT")
    label, raw_path = value.split("=", 1)
    root = Path(raw_path).expanduser().resolve()
    summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
    if summary.get("status") != "DONE":
        raise RuntimeError(f"Run {label} is not DONE")
    return {
        "label": label,
        "output_root": str(root),
        "contract_id": summary["contract_id"],
        "total": int(summary["total"]),
        "successful": int(summary["successful"]),
        "failed": int(summary["failed"]),
        "wall_time_seconds": float(summary["wall_time_seconds"]),
        "molecules_per_second": float(summary["molecules_per_second"]),
        "stage_timing_seconds": summary["stage_timing_seconds"],
        "operational": summary["operational"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    rows = [load_run(value) for value in args.run]
    contracts = {row["contract_id"] for row in rows}
    if len(contracts) != 1:
        raise RuntimeError("Cannot rank different scientific contracts")
    eligible = [row for row in rows if row["failed"] == 0 and row["successful"] == row["total"]]
    ranking = sorted(eligible, key=lambda row: (-row["molecules_per_second"], row["label"]))
    for rank, row in enumerate(ranking, 1):
        row["throughput_rank"] = rank
    output = args.output_root.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    document = {"schema_version": 1, "generated_at": utc_now(), "runs": rows, "ranking": [row["label"] for row in ranking]}
    atomic_write_json(output / "performance_summary.json", document)
    columns = ["label", "total", "successful", "failed", "wall_time_seconds", "molecules_per_second", "throughput_rank"]
    csv_path = output / "performance_summary.csv"
    temporary = csv_path.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})
    temporary.replace(csv_path)
    lines = ["# Demiurge performance summary", "", "| Rank | Run | Wall (s) | molecules/s | Failed |", "|---:|---|---:|---:|---:|"]
    for row in sorted(rows, key=lambda item: item.get("throughput_rank", 10**9)):
        lines.append(f"| {row.get('throughput_rank', '-')} | {row['label']} | {row['wall_time_seconds']:.3f} | {row['molecules_per_second']:.4f} | {row['failed']} |")
    if ranking:
        lines.extend(["", f"Recommended fastest exact/QC-clean run: **{ranking[0]['label']}**."])
    else:
        lines.extend(["", "No run is recommendation-eligible because scientific/QC failures were present."])
    atomic_write_text(output / "performance_summary.md", "\n".join(lines) + "\n")
    print(json.dumps(document, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
