"""Shared NMR V2 scientific pipeline used by local and SLURM backends."""

from __future__ import annotations

import hashlib
import json
import multiprocessing
import os
import shutil
import signal
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from . import predictor
from .computation_cache import ComputationCache
from .bucketing import (
    bucket_shifts,
    parse_prediction_csv,
    validate_prediction_records,
)
from .contracts import (
    C_BINS,
    C_MAX,
    C_MIN,
    ECFP_BITS,
    ECFP_RADIUS,
    ECFP_USE_CHIRALITY,
    H_BINS,
    H_MAX,
    H_MIN,
    MODE_FEATURE_DIMENSIONS,
    MOLECULAR_IDENTITY_CONTRACT,
    NMR_REPRESENTATION_VERSION,
    RECORD_ID_CONTRACT,
    object_sha256,
    scientific_contract,
    verify_predictor_artifacts,
)
from .io_utils import cleanup_owned_scratch, create_owned_scratch
from .io import create_input_reader, create_output_writer
from .io.base import ColumnSelector
from .java_heap import DEFAULT_JAVA_HEAP, normalize_java_heap
from .preparation import (
    MolecularIdentity,
    PreparationResult,
    _worker_init,
    molecular_identity_v2,
    prepare_batch,
)
from .run_state import (
    atomic_write_json,
    file_sha256,
    initial_checkpoint,
    input_identity,
    read_json,
    utc_now,
    validate_resume,
    write_progress,
)


_STOP_SIGNAL: int | None = None


@dataclass(frozen=True)
class RunConfig:
    input_path: Path
    mode: str
    output_root: Path
    temp_root: Path
    label_column: ColumnSelector = 3
    prep_workers: int = 4
    java_threads: int = 2
    java_heap: str = DEFAULT_JAVA_HEAP
    batch_size: int = 500
    java_lifecycle: str = "persistent"
    max_attempts: int = 3
    retain_scientific_artifacts: bool = False
    backend: str = "local"
    resume: bool = False
    canonical_input_path: Path | None = None
    input_format: str = "csv"
    input_table: str | None = None
    input_query: str | None = None
    id_column: str = "MOLECULE_NAME"
    smiles_column: str = "SMILES"
    output_format: str = "csv"
    output_db: Path | None = None
    output_table: str = "demiurge_features"
    metadata_table: str = "demiurge_metadata"
    overwrite_output: bool = False
    include_murcko: bool = False

    def validated(self) -> "RunConfig":
        mode_aliases = {"1h": "1H", "13c": "13C", "fp": "FP", "hybrid": "hybrid", "total": "total"}
        mode = mode_aliases.get(str(self.mode).lower())
        if mode is None:
            raise ValueError("mode must be one of 1H, 13C, FP, hybrid or total")
        for name in ("prep_workers", "java_threads", "batch_size", "max_attempts"):
            value = getattr(self, name)
            if isinstance(value, bool) or int(value) <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if isinstance(self.label_column, int) and self.label_column <= 0:
            raise ValueError("label_column position must be positive")
        if not isinstance(self.label_column, (int, str)) or str(self.label_column).strip() == "":
            raise ValueError("label_column must be a one-based position or column name")
        if self.java_lifecycle not in {"per-batch", "persistent"}:
            raise ValueError("java_lifecycle must be per-batch or persistent")
        if self.backend not in {"local", "slurm-worker"}:
            raise ValueError("backend must be local or slurm-worker")
        if self.input_format not in {"csv", "sqlite"}:
            raise ValueError("input_format must be csv or sqlite")
        if self.output_format not in {"csv", "sqlite"}:
            raise ValueError("output_format must be csv or sqlite")
        input_path = self.input_path.expanduser().resolve()
        if not input_path.is_file():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")
        output_db = self.output_db.expanduser().resolve() if self.output_db is not None else None
        if output_db is not None and output_db == input_path:
            raise ValueError("SQLite output_db must not be the input database")
        return replace(
            self,
            input_path=input_path,
            mode=mode,
            output_root=self.output_root.expanduser().resolve(),
            temp_root=self.temp_root.expanduser().resolve(),
            java_heap=normalize_java_heap(self.java_heap),
            canonical_input_path=(
                self.canonical_input_path.expanduser().resolve()
                if self.canonical_input_path is not None
                else input_path
            ),
            output_db=output_db,
        )


def install_signal_handlers() -> None:
    def handler(signum: int, _frame: Any) -> None:
        global _STOP_SIGNAL
        _STOP_SIGNAL = signum
        predictor.shutdown_persistent_java_processors(f"python-signal-{signum}")

    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)


def _scientific_config(config: RunConfig, label_name: str) -> dict[str, Any]:
    from rdkit import rdBase

    contract = scientific_contract(config.mode)
    artifacts = (
        verify_predictor_artifacts(Path(__file__).resolve().parent.parent)
        if config.mode != "FP"
        else {}
    )
    return {
        "contract": contract,
        "contract_sha256": object_sha256(contract),
        "mode": config.mode,
        "label_column": config.label_column,
        "label_name": label_name,
        "rdkit_version": rdBase.rdkitVersion,
        "predictor_artifact_sha256": artifacts,
        "molecular_identity_contract": MOLECULAR_IDENTITY_CONTRACT,
    }


def _create_ecfp_generator():
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

    return GetMorganGenerator(
        radius=ECFP_RADIUS,
        fpSize=ECFP_BITS,
        includeChirality=ECFP_USE_CHIRALITY,
    )


def _ecfp(smiles: str, generator: Any) -> list[int]:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit rejected SMILES for ECFP4")
    return [int(value) for value in generator.GetFingerprint(mol).ToBitString()]


def _nucleus_vector(path: Path, nucleus: str) -> tuple[list[int], dict[str, int]]:
    records = parse_prediction_csv(path, nucleus)
    shifts = validate_prediction_records(records, nucleus)
    if not shifts:
        raise ValueError(f"Empty raw {nucleus} prediction")
    if nucleus == "1H":
        vector, inside, outside = bucket_shifts(shifts, H_MIN, H_MAX, H_BINS)
    else:
        vector, inside, outside = bucket_shifts(shifts, C_MIN, C_MAX, C_BINS)
    if int(np.count_nonzero(vector)) == 0:
        raise ValueError(f"All-zero {nucleus} NMR vector")
    return [int(value) for value in vector.tolist()], {
        "total": len(shifts), "in_range": inside, "out_of_range": outside,
    }


def _failure(record: dict[str, Any], stage: str, error: Exception | str) -> dict[str, Any]:
    message = str(error)
    error_type = type(error).__name__ if isinstance(error, Exception) else "InputQCError"
    return {
        "source_index": record["source_index"],
        "record_id": record["record_id"],
        "internal_id": record["internal_id"],
        "molecule_name": record["molecule_name"],
        "smiles": record["smiles"],
        "label": record["label"],
        "status": "FAILED",
        "failure_stage": stage,
        "failure_type": error_type,
        "failure_message": message,
    }


def _process_batch(
    records: list[dict[str, Any]],
    config: RunConfig,
    scratch_batch: Path,
    prep_pool: Any,
    feature_cache: ComputationCache,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, float]]:
    timings = {"preparation": 0.0, "java_1h": 0.0, "java_13c": 0.0, "features": 0.0}
    mol_dir = scratch_batch / "mols"
    raw_h = scratch_batch / "raw_1h"
    raw_c = scratch_batch / "raw_13c"
    for directory in (mol_dir, raw_h, raw_c):
        directory.mkdir(parents=True, exist_ok=True)

    metadata: list[dict[str, Any]] = []
    good_input = []
    for record in records:
        if record["input_error"]:
            metadata.append(_failure(record, "INPUT_QC", record["input_error"]))
        else:
            good_input.append(record)

    needs_nmr = config.mode != "FP"
    preparation_by_id: dict[str, PreparationResult] = {}
    identity_by_id: dict[str, MolecularIdentity] = {}
    if needs_nmr and good_input:
        started = time.perf_counter()
        prepared = prepare_batch(good_input, mol_dir, pool=prep_pool)
        timings["preparation"] = time.perf_counter() - started
        preparation_by_id = {item.internal_id: item for item in prepared}
        for record in good_input:
            result = preparation_by_id[record["internal_id"]]
            if not result.successful:
                metadata.append(_failure(
                    record,
                    "PREPARATION",
                    f"{result.error_type}: {result.error_message}",
                ))
            else:
                identity_by_id[record["internal_id"]] = MolecularIdentity(
                    canonical_smiles=str(result.canonical_smiles),
                    identity_smiles=str(result.identity_smiles),
                    identity_sha256=str(result.identity_sha256),
                    murcko_smiles=str(result.murcko_smiles),
                    murcko_id=str(result.murcko_id),
                )
    elif good_input:
        for record in good_input:
            try:
                identity_by_id[record["internal_id"]] = molecular_identity_v2(record["smiles"])
            except Exception as exc:
                metadata.append(_failure(record, "PREPARATION", exc))

    prepared_records = [
        record for record in good_input
        if record["internal_id"] in identity_by_id
    ]
    groups: dict[str, list[dict[str, Any]]] = {}
    for record in prepared_records:
        key = identity_by_id[record["internal_id"]].identity_sha256
        groups.setdefault(key, []).append(record)

    uncached_representatives: list[dict[str, Any]] = []
    cached_by_key: dict[str, tuple[list[int], dict[str, Any]]] = {}
    for key, group in groups.items():
        identity = identity_by_id[group[0]["internal_id"]]
        cached = None if config.retain_scientific_artifacts else feature_cache.get(
            key, identity.identity_smiles
        )
        if cached is None:
            uncached_representatives.append(group[0])
        else:
            cached_by_key[key] = (cached.vector, cached.diagnostics)

    representative_ids = {item["internal_id"] for item in uncached_representatives}
    if needs_nmr:
        for record in prepared_records:
            if record["internal_id"] not in representative_ids:
                Path(str(preparation_by_id[record["internal_id"]].mol_path)).unlink(missing_ok=True)

    if needs_nmr and uncached_representatives:
        if config.mode in {"1H", "hybrid", "total"}:
            started = time.perf_counter()
            result = predictor.run_java_batch_processor(
                mol_dir,
                "1H",
                java_threads=config.java_threads,
                java_heap=config.java_heap,
                output_directory=raw_h,
            )
            timings["java_1h"] = time.perf_counter() - started
            if result is None:
                raise RuntimeError("1H Java batch subprocess failed")
        if config.mode in {"13C", "hybrid", "total"}:
            started = time.perf_counter()
            result = predictor.run_java_batch_processor(
                mol_dir,
                "13C",
                java_threads=config.java_threads,
                java_heap=config.java_heap,
                output_directory=raw_c,
            )
            timings["java_13c"] = time.perf_counter() - started
            if result is None:
                raise RuntimeError("13C Java batch subprocess failed")

    generator = _create_ecfp_generator() if config.mode in {"FP", "total"} else None
    feature_rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    computed_by_key = dict(cached_by_key)
    error_by_key: dict[str, Exception] = {}
    for record in uncached_representatives:
        identity = identity_by_id[record["internal_id"]]
        key = identity.identity_sha256
        try:
            h_vector: list[int] = []
            c_vector: list[int] = []
            diagnostics: dict[str, Any] = {}
            if config.mode in {"1H", "hybrid", "total"}:
                path = raw_h / f"{record['internal_id']}.csv"
                if not path.is_file():
                    raise FileNotFoundError(f"Missing 1H predictor CSV {path}")
                h_vector, diagnostics["1H"] = _nucleus_vector(path, "1H")
            if config.mode in {"13C", "hybrid", "total"}:
                path = raw_c / f"{record['internal_id']}.csv"
                if not path.is_file():
                    raise FileNotFoundError(f"Missing 13C predictor CSV {path}")
                c_vector, diagnostics["13C"] = _nucleus_vector(path, "13C")
            fp_vector = _ecfp(record["smiles"], generator) if generator is not None else []
            if config.mode == "1H":
                vector = h_vector
            elif config.mode == "13C":
                vector = c_vector
            elif config.mode == "hybrid":
                vector = h_vector + c_vector
            elif config.mode == "FP":
                vector = fp_vector
            else:
                vector = h_vector + c_vector + fp_vector
            if len(vector) != MODE_FEATURE_DIMENSIONS[config.mode]:
                raise RuntimeError(f"Feature dimension mismatch: {len(vector)}")
            computed_by_key[key] = (vector, diagnostics)
            feature_cache.put(key, identity.identity_smiles, vector, diagnostics)
        except Exception as exc:
            error_by_key[key] = exc

    for key, group in groups.items():
        identity = identity_by_id[group[0]["internal_id"]]
        if key in error_by_key:
            for record in group:
                metadata.append(_failure(record, "FEATURE_ASSEMBLY", error_by_key[key]))
            continue
        vector, diagnostics = computed_by_key[key]
        representative = group[0]
        if needs_nmr and len(group) > 1:
            representative_mol = Path(str(preparation_by_id[representative["internal_id"]].mol_path))
            for record in group[1:]:
                target_mol = mol_dir / f"{record['internal_id']}.mol"
                if not target_mol.exists() and representative_mol.exists():
                    shutil.copyfile(representative_mol, target_mol)
                for source_dir in (raw_h, raw_c):
                    source_csv = source_dir / f"{representative['internal_id']}.csv"
                    target_csv = source_dir / f"{record['internal_id']}.csv"
                    if source_csv.exists() and not target_csv.exists():
                        shutil.copyfile(source_csv, target_csv)
        for record in group:
            feature_rows.append({
                "source_index": record["source_index"],
                "record_id": record["record_id"],
                "molecule_name": record["molecule_name"],
                "label": record["label"],
                "murcko_smiles": identity.murcko_smiles,
                "murcko_id": identity.murcko_id,
                "vector": vector,
            })
            metadata.append({
                "source_index": record["source_index"],
                "record_id": record["record_id"],
                "internal_id": record["internal_id"],
                "molecule_name": record["molecule_name"],
                "smiles": record["smiles"],
                "canonical_smiles": identity.canonical_smiles,
                "molecular_identity_smiles": identity.identity_smiles,
                "molecular_identity_sha256": identity.identity_sha256,
                "murcko_smiles": identity.murcko_smiles if config.include_murcko else None,
                "murcko_id": identity.murcko_id if config.include_murcko else None,
                "status": "SUCCESS",
                "nmr_diagnostics": diagnostics,
                "feature_sha256": hashlib.sha256(
                    json.dumps(vector, separators=(",", ":")).encode("ascii")
                ).hexdigest(),
            })
            feature_cache.observe_record(
                record_id=record["record_id"],
                molecule_name=record["molecule_name"],
                raw_smiles=record["smiles"],
                label=record["label"],
                identity_sha256=identity.identity_sha256,
            )
    timings["features"] = time.perf_counter() - started
    feature_rows.sort(key=lambda item: int(item["source_index"]))
    metadata.sort(key=lambda item: int(item["source_index"]))
    return feature_rows, metadata, timings


def _check_stop() -> None:
    if _STOP_SIGNAL is not None:
        raise InterruptedError(f"Demiurge interrupted by signal {_STOP_SIGNAL}")


def run_pipeline(config: RunConfig) -> dict[str, Any]:
    global _STOP_SIGNAL
    _STOP_SIGNAL = None
    config = config.validated()
    install_signal_handlers()
    config.output_root.mkdir(parents=True, exist_ok=True)
    reader = create_input_reader(
        config.input_path,
        input_format=config.input_format,
        input_table=config.input_table,
        input_query=config.input_query,
        id_column=config.id_column,
        smiles_column=config.smiles_column,
        label_column=config.label_column,
    )
    description = reader.describe()
    scientific = _scientific_config(config, description.label_name)
    writer = create_output_writer(
        output_format=config.output_format,
        output_root=config.output_root,
        output_db=config.output_db,
        output_table=config.output_table,
        metadata_table=config.metadata_table,
        input_stem=config.input_path.stem,
        mode=config.mode,
        feature_contract=scientific["contract"],
        overwrite_output=config.overwrite_output,
        include_murcko=config.include_murcko,
    )
    io_config = {"input": description.configuration, "output": writer.configuration()}
    identity = input_identity(config.input_path)
    run_id = hashlib.sha256(
        (identity["sha256"] + object_sha256(scientific)).encode("ascii")
    ).hexdigest()[:16]
    manifest_path = config.output_root / "run_manifest.json"
    checkpoint_path = config.output_root / "checkpoint.json"
    if manifest_path.exists() and not config.resume:
        raise RuntimeError(f"Output root already contains a run manifest: {manifest_path}")

    if config.resume:
        if not manifest_path.is_file() or not checkpoint_path.is_file():
            raise RuntimeError("Resume requires run_manifest.json and checkpoint.json")
        manifest = read_json(manifest_path)
        checkpoint = read_json(checkpoint_path)
        validate_resume(
            checkpoint,
            input_info=identity,
            scientific_config=scientific,
            io_config=io_config,
        )
        checkpoint["attempt"] = int(checkpoint.get("attempt", 1)) + 1
        if checkpoint["attempt"] > int(config.max_attempts):
            raise RuntimeError("Maximum run attempts exhausted")
        writer.initialize(resume=True)
    else:
        manifest = {
            "manifest_schema_version": 1,
            "run_id": run_id,
            "input_path": str(config.canonical_input_path or config.input_path),
            "input_identity": identity,
            "mode": config.mode,
            "label_column": config.label_column,
            "label_name": description.label_name,
            "scientific_config": scientific,
            "scientific_config_sha256": object_sha256(scientific),
            "io_config": io_config,
            "record_id_contract": RECORD_ID_CONTRACT,
            "operational_initial": {
                "prep_workers": config.prep_workers,
                "java_threads": config.java_threads,
                "java_heap": config.java_heap,
                "batch_size": config.batch_size,
                "java_lifecycle": config.java_lifecycle,
                "max_attempts": config.max_attempts,
                "retain_scientific_artifacts": config.retain_scientific_artifacts,
                "include_murcko": config.include_murcko,
            },
            "created_at": utc_now(),
        }
        writer.initialize(resume=False)
        atomic_write_json(manifest_path, manifest)
        checkpoint = initial_checkpoint(
            run_id=run_id,
            input_info=identity,
            scientific_config=scientific,
            total_rows=description.total_rows,
            io_config=io_config,
        )
    checkpoint.update({"status": "RUNNING", "heartbeat": utc_now(), "pid": os.getpid(), "backend": config.backend})
    atomic_write_json(checkpoint_path, checkpoint)
    write_progress(config.output_root, checkpoint)

    os.environ[predictor.JAVA_LIFECYCLE_ENV] = config.java_lifecycle
    os.environ[predictor.JAVA_PREDICTOR_MODE_ENV] = predictor.PREDICTOR_MODE_THREAD_LOCAL
    os.environ["SPECTRAPRINTS_UNIFIED_PROFILE"] = "1"
    os.environ[predictor.JAVA_DIAGNOSTICS_DIR_ENV] = str(config.output_root / "diagnostics" / "java")
    scratch, owner_token = create_owned_scratch(config.temp_root, run_id)
    feature_cache = ComputationCache(scratch / "molecular_feature_cache.sqlite")
    prep_pool = None
    batch_iterator = None
    total_started = time.perf_counter()
    aggregate_timings = {"preparation": 0.0, "java_1h": 0.0, "java_13c": 0.0, "features": 0.0}
    try:
        if config.mode != "FP" and config.prep_workers > 1:
            prep_pool = multiprocessing.get_context("spawn").Pool(
                processes=config.prep_workers,
                initializer=_worker_init,
            )
        start_index = int(checkpoint.get("next_row_index", 0))
        batch_index = int(checkpoint.get("committed_batches", 0))
        batch_iterator = reader.iter_batches(config.batch_size, start_index)
        for offset, selected in batch_iterator:
            _check_stop()
            for record in selected:
                row_identity = json.dumps(
                    [record["molecule_name"], record["smiles"], record["label"]],
                    ensure_ascii=False, allow_nan=False, separators=(",", ":"),
                ).encode("utf-8")
                record["record_id"] = (
                    f"{hashlib.sha256(row_identity).hexdigest()[:16]}-"
                    f"R{int(record['source_index']) + 1:012d}"
                )
            last_error: Exception | None = None
            for batch_attempt in range(1, config.max_attempts + 1):
                scratch_batch = scratch / f"batch_{batch_index:08d}_try_{batch_attempt:02d}"
                if scratch_batch.exists():
                    shutil.rmtree(scratch_batch)
                scratch_batch.mkdir()
                try:
                    feature_rows, metadata, timings = _process_batch(
                        selected, config, scratch_batch, prep_pool, feature_cache
                    )
                    writer.commit_batch(
                        batch_index,
                        feature_rows,
                        metadata,
                        timings,
                        scratch_batch,
                        config.retain_scientific_artifacts,
                    )
                    last_error = None
                    break
                except (OSError, RuntimeError, TimeoutError) as exc:
                    last_error = exc
                    predictor.shutdown_persistent_java_processors("transient-batch-failure")
                    if batch_attempt >= config.max_attempts:
                        raise
                finally:
                    # A persistent JVM may retain native/CDK handles to the
                    # current MOL files until shutdown (observable on Windows).
                    # Keep such batches inside the owned run scratch and remove
                    # the whole marker-owned tree after JVM shutdown below.
                    if scratch_batch.exists() and config.java_lifecycle == "per-batch":
                        shutil.rmtree(scratch_batch)
            if last_error is not None:
                raise last_error
            for key, value in timings.items():
                aggregate_timings[key] += float(value)
            successful = sum(1 for item in metadata if item["status"] == "SUCCESS")
            failed = sum(1 for item in metadata if item["status"] == "FAILED")
            batch_index += 1
            checkpoint.update({
                "status": "RUNNING",
                "next_row_index": offset + len(selected),
                "committed_batches": batch_index,
                "successful": int(checkpoint["successful"]) + successful,
                "failed": int(checkpoint["failed"]) + failed,
                "heartbeat": utc_now(),
                "last_batch_timing_seconds": timings,
            })
            atomic_write_json(checkpoint_path, checkpoint)
            write_progress(config.output_root, checkpoint)

        final, failures = writer.finalize(
            total=description.total_rows,
            successful=int(checkpoint["successful"]),
            failed=int(checkpoint["failed"]),
        )
        wall = time.perf_counter() - total_started
        checkpoint.update({"status": "DONE", "heartbeat": utc_now(), "finished_at": utc_now(), "runtime_seconds": wall})
        atomic_write_json(checkpoint_path, checkpoint)
        summary = {
            "summary_schema_version": 1,
            "status": "DONE",
            "run_id": run_id,
            "contract_id": scientific["contract"]["contract_id"],
            "nmr_representation_version": NMR_REPRESENTATION_VERSION if config.mode != "FP" else None,
            "input_identity": identity,
            "total": description.total_rows,
            "successful": checkpoint["successful"],
            "failed": checkpoint["failed"],
            "wall_time_seconds": wall,
            "molecules_per_second": (description.total_rows / wall if wall else None),
            "stage_timing_seconds": aggregate_timings,
            "record_audit": feature_cache.audit(),
            "final_output": str(final),
            "final_output_sha256": file_sha256(final),
            "failures": str(failures),
            "io": io_config,
            "operational": {
                "backend": config.backend,
                "batch_size": config.batch_size,
                "prep_workers": config.prep_workers,
                "java_threads": config.java_threads,
                "java_heap": config.java_heap,
                "java_lifecycle": config.java_lifecycle,
                "runtime_temp_root": str(scratch),
            },
            "completed_at": utc_now(),
        }
        atomic_write_json(config.output_root / "summary.json", summary)
        write_progress(config.output_root, checkpoint)
        return summary
    except InterruptedError as exc:
        checkpoint.update({"status": "INTERRUPTED", "heartbeat": utc_now(), "failure_type": type(exc).__name__, "failure_message": str(exc)})
        atomic_write_json(checkpoint_path, checkpoint)
        write_progress(config.output_root, checkpoint)
        raise
    except Exception as exc:
        checkpoint.update({"status": "FAILED", "heartbeat": utc_now(), "failure_type": type(exc).__name__, "failure_message": str(exc)})
        atomic_write_json(checkpoint_path, checkpoint)
        write_progress(config.output_root, checkpoint)
        raise
    finally:
        close_iterator = getattr(batch_iterator, "close", None)
        if callable(close_iterator):
            close_iterator()
        if prep_pool is not None:
            prep_pool.terminate()
            prep_pool.join()
        predictor.shutdown_persistent_java_processors("pipeline-finally")
        feature_cache.close()
        cleanup_owned_scratch(config.temp_root, scratch, owner_token)


def resume_pipeline(output_root: Path, temp_root: Path | None = None, **overrides: Any) -> dict[str, Any]:
    root = output_root.expanduser().resolve()
    manifest = read_json(root / "run_manifest.json")
    initial = manifest.get("operational_initial") or {}
    io_config = manifest.get("io_config") or {
        "input": {
            "input_format": "csv",
            "id_column": "MOLECULE_NAME",
            "smiles_column": "SMILES",
            "label_column": manifest["label_column"],
            "label_name": manifest["label_name"],
        },
        "output": {"output_format": "csv"},
    }
    input_config = io_config["input"]
    output_config = io_config["output"]
    values = {
        "input_path": Path(overrides.get("input_path") or manifest["input_path"]),
        "mode": manifest["mode"],
        "output_root": root,
        "temp_root": temp_root or Path(initial.get("temp_root") or root / "tmp"),
        "label_column": input_config["label_column"],
        "prep_workers": int(overrides.get("prep_workers") or initial.get("prep_workers", 4)),
        "java_threads": int(overrides.get("java_threads") or initial.get("java_threads", 2)),
        "java_heap": str(overrides.get("java_heap") or initial.get("java_heap", DEFAULT_JAVA_HEAP)),
        "batch_size": int(overrides.get("batch_size") or initial.get("batch_size", 500)),
        "java_lifecycle": str(overrides.get("java_lifecycle") or initial.get("java_lifecycle", "persistent")),
        "max_attempts": int(initial.get("max_attempts", 3)),
        "retain_scientific_artifacts": bool(initial.get("retain_scientific_artifacts", False)),
        "backend": str(overrides.get("backend") or "local"),
        "resume": True,
        "canonical_input_path": Path(manifest["input_path"]),
        "input_format": input_config.get("input_format", "csv"),
        "input_table": input_config.get("input_table"),
        "input_query": input_config.get("input_query"),
        "id_column": input_config.get("id_column", "MOLECULE_NAME"),
        "smiles_column": input_config.get("smiles_column", "SMILES"),
        "output_format": output_config.get("output_format", "csv"),
        "output_db": Path(output_config["output_db"]) if output_config.get("output_db") else None,
        "output_table": output_config.get("output_table", "demiurge_features"),
        "metadata_table": output_config.get("metadata_table", "demiurge_metadata"),
        "include_murcko": bool(initial.get("include_murcko", False)),
    }
    return run_pipeline(RunConfig(**values))


def read_status(output_root: Path) -> tuple[str, Path]:
    root = output_root.expanduser().resolve()
    checkpoint = read_json(root / "checkpoint.json")
    path = write_progress(root, checkpoint)
    return path.read_text(encoding="utf-8"), path
