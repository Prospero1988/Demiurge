# Project invariants

- `SPECTRAPRINTS_NMR_V2` is the active production NMR representation. Do not reintroduce ETKDG, CoordGen or OpenBabel into the active path.
- The mandatory deployment reference is the hash-pinned frozen NMR V2 fixture in `validation/frozen_nmr_v2_expected`; normal validation must remain standalone. `screen_SPECTRAprints` is the historical source and an optional development-only cross-repository gate when both checkouts are available.
- Feature orders are frozen: 1H=200; 13C=200; hybrid=H|C=400; FP=ECFP4 2048; total=H|C|ECFP4=2448.
- Operational settings (backend, paths, staging, batching, worker counts, Java thread count/heap/lifecycle) must not enter scientific identity or alter output.
- Local and SLURM execution must call the same scientific core. Do not duplicate scientific logic in orchestration.
- Predictor JAR hashes are part of the scientific contract. Never silently accept a changed artifact.
- Preserve fail-closed checkpoint compatibility, marker-owned scratch cleanup, finite retries and atomic batch commits.
- After every accepted functional, architectural, orchestration or performance change, update `README.md`. For performance work, preserve historical benchmark results rather than replacing them with only the newest result.
