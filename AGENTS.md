# Project invariants

- `SPECTRAPRINTS_NMR_V2` is the active production NMR representation. Do not reintroduce ETKDG, CoordGen or OpenBabel into the active path.
- The exact scientific reference is the validated implementation in `screen_SPECTRAprints`. Any preparation, MOL-byte, raw-spectrum or bucket change requires the cross-repository zero-tolerance parity gate.
- Feature orders are frozen: 1H=200; 13C=200; hybrid=H|C=400; FP=ECFP4 2048; total=H|C|ECFP4=2448.
- Operational settings (backend, paths, staging, batching, worker counts, Java thread count/heap/lifecycle) must not enter scientific identity or alter output.
- Local and SLURM execution must call the same scientific core. Do not duplicate scientific logic in orchestration.
- Predictor JAR hashes are part of the scientific contract. Never silently accept a changed artifact.
- Preserve fail-closed checkpoint compatibility, marker-owned scratch cleanup, finite retries and atomic batch commits.
- After every accepted functional, architectural, orchestration or performance change, update `README.md`. For performance work, preserve historical benchmark results rather than replacing them with only the newest result.
