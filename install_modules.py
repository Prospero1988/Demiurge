#!/usr/bin/env python3
"""Fail-closed dependency preflight for the production Demiurge runtime."""

from __future__ import annotations

import importlib.util
import json
import sys

from demiurge_bin.preflight import run_preflight


REQUIRED_MODULES = ("numpy", "pandas", "rdkit")


def main() -> int:
    missing_modules = [name for name in REQUIRED_MODULES if importlib.util.find_spec(name) is None]
    if missing_modules:
        if missing_modules:
            print("Missing Python modules: " + ", ".join(missing_modules), file=sys.stderr)
        print("Create/update the Conda environment from conda_environment.yml.", file=sys.stderr)
        return 1
    try:
        report = run_preflight()
    except Exception as error:
        print(f"Demiurge production dependency preflight: FAIL: {error}", file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2, sort_keys=True))
    print("Demiurge production dependency preflight: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
