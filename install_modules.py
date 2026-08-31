#!/usr/bin/env python3
"""Fail-closed dependency preflight for the production Demiurge runtime."""

from __future__ import annotations

import importlib.util
import shutil
import sys


REQUIRED_MODULES = ("numpy", "pandas", "rdkit")
REQUIRED_EXECUTABLES = ("java", "javac")


def main() -> int:
    missing_modules = [name for name in REQUIRED_MODULES if importlib.util.find_spec(name) is None]
    missing_commands = [name for name in REQUIRED_EXECUTABLES if shutil.which(name) is None]
    if missing_modules or missing_commands:
        if missing_modules:
            print("Missing Python modules: " + ", ".join(missing_modules), file=sys.stderr)
        if missing_commands:
            print("Missing commands: " + ", ".join(missing_commands), file=sys.stderr)
        print("Create conda_environment.yml and install a JDK; OpenBabel is not required.", file=sys.stderr)
        return 1
    print("Demiurge production dependency preflight: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
