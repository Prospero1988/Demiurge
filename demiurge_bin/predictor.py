import hashlib
import atexit
import base64
import itertools
import json
import os
import queue
import subprocess
import platform
import threading
import time
import tempfile
from pathlib import Path

from .java_heap import DEFAULT_JAVA_HEAP, normalize_java_heap
from .contracts import verify_predictor_artifacts


_LAST_JAVA_COMMANDS: dict[str, tuple[str, ...]] = {}
_DIAGNOSTIC_SEQUENCE = itertools.count(1)

JAVA_DIAGNOSTICS_DIR_ENV = "SPECTRAPRINTS_JAVA_DIAGNOSTICS_DIR"
JAVA_JFR_ENV = "SPECTRAPRINTS_JAVA_JFR"
JAVA_PREDICTOR_MODE_ENV = "SPECTRAPRINTS_JAVA_PREDICTOR_MODE"
JAVA_LIFECYCLE_ENV = "SPECTRAPRINTS_JAVA_LIFECYCLE"
JAVA_PROTOCOL_TIMEOUT_ENV = "SPECTRAPRINTS_JAVA_PROTOCOL_TIMEOUT_SECONDS"
JAVA_BUILD_DIR_ENV = "DEMIURGE_JAVA_BUILD_DIR"
PREDICTOR_MODE_THREAD_LOCAL = "thread-local"
PREDICTOR_MODE_PER_MOLECULE = "per-molecule"
JAVA_LIFECYCLE_PER_BATCH = "per-batch"
JAVA_LIFECYCLE_PERSISTENT = "persistent"
DEFAULT_JAVA_LIFECYCLE = JAVA_LIFECYCLE_PERSISTENT
_VALID_JAVA_LIFECYCLES = {JAVA_LIFECYCLE_PER_BATCH, JAVA_LIFECYCLE_PERSISTENT}
_PROTOCOL_PREFIX = "SPECTRAPRINTS_C3"
_PROTOCOL_VERSION = "1"
_PERSISTENT_PROCESSORS: dict[tuple[str, int, str, str], "_PersistentJavaProcessor"] = {}
_PERSISTENT_LOCK = threading.RLock()
_VALID_PREDICTOR_MODES = {
    PREDICTOR_MODE_THREAD_LOCAL,
    PREDICTOR_MODE_PER_MOLECULE,
}


def _environment_flag(name: str) -> bool:
    value = os.environ.get(name, "").strip().lower()
    if value in ("", "0", "false", "no", "off"):
        return False
    if value in ("1", "true", "yes", "on"):
        return True
    raise ValueError(f"{name} must be one of 0/1, false/true, no/yes, or off/on")


def _predictor_mode_from_environment() -> str:
    value = os.environ.get(
        JAVA_PREDICTOR_MODE_ENV, PREDICTOR_MODE_THREAD_LOCAL
    ).strip().lower()
    if value not in _VALID_PREDICTOR_MODES:
        choices = ", ".join(sorted(_VALID_PREDICTOR_MODES))
        raise ValueError(f"{JAVA_PREDICTOR_MODE_ENV} must be one of: {choices}")
    return value


def java_lifecycle_from_environment() -> str:
    value = os.environ.get(
        JAVA_LIFECYCLE_ENV, DEFAULT_JAVA_LIFECYCLE
    ).strip().lower()
    if value not in _VALID_JAVA_LIFECYCLES:
        choices = ", ".join(sorted(_VALID_JAVA_LIFECYCLES))
        raise ValueError(f"{JAVA_LIFECYCLE_ENV} must be one of: {choices}")
    return value


def _protocol_timeout_from_environment() -> float:
    raw = os.environ.get(JAVA_PROTOCOL_TIMEOUT_ENV, "7200").strip()
    try:
        value = float(raw)
    except ValueError as error:
        raise ValueError(f"{JAVA_PROTOCOL_TIMEOUT_ENV} must be a positive number") from error
    if not value > 0.0:
        raise ValueError(f"{JAVA_PROTOCOL_TIMEOUT_ENV} must be a positive number")
    return value


def _diagnostic_paths(predictor: str) -> dict[str, Path] | None:
    raw_directory = os.environ.get(JAVA_DIAGNOSTICS_DIR_ENV, "").strip()
    jfr_enabled = _environment_flag(JAVA_JFR_ENV)
    if not raw_directory:
        if jfr_enabled:
            raise ValueError(
                f"{JAVA_JFR_ENV}=1 requires {JAVA_DIAGNOSTICS_DIR_ENV}"
            )
        return None

    directory = Path(raw_directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    sequence = next(_DIAGNOSTIC_SEQUENCE)
    scope_name = (
        os.environ.get("SLURM_JOB_ID")
        or os.environ.get("HOSTNAME")
        or os.environ.get("COMPUTERNAME")
        or "host"
    )
    scope = f"{scope_name}_p{os.getpid()}"
    safe_scope = "".join(
        character if character.isalnum() or character in ("-", "_") else "_"
        for character in scope
    )
    stem = f"java_{predictor.lower()}_{safe_scope}_{sequence:06d}"
    paths = {
        "directory": directory,
        "java": directory / f"{stem}.java.json",
        "launcher": directory / f"{stem}.launcher.json",
        "gc": directory / f"{stem}.gc.log",
        "resources": directory / f"{stem}.resources.txt",
    }
    if jfr_enabled:
        paths["jfr"] = directory / f"{stem}.jfr"
    return paths


def _write_launcher_diagnostics(
    path: Path,
    *,
    predictor: str,
    predictor_mode: str,
    java_threads: int,
    java_heap: str,
    command: list[str],
    elapsed_seconds: float,
    returncode: int | None,
    java_metrics_path: Path,
    lifecycle: str = JAVA_LIFECYCLE_PER_BATCH,
    startup_seconds: float | None = None,
    batches: list[dict[str, object]] | None = None,
    restarts: int = 0,
    protocol_failures: int = 0,
    protocol_timeouts: int = 0,
    shutdown_reason: str | None = None,
    jvm_startup_count: int = 1,
    startup_measurement: str | None = None,
) -> None:
    document = {
        "schema_version": 1,
        "predictor": predictor,
        "predictor_mode": predictor_mode,
        "java_threads": java_threads,
        "java_heap": java_heap,
        "subprocess_wall_seconds": elapsed_seconds,
        "returncode": returncode,
        "java_metrics_path": str(java_metrics_path),
        "command": command,
        "java_lifecycle": lifecycle,
        "jvm_startup_count": jvm_startup_count,
        "jvm_startup_wall_seconds": startup_seconds,
        "jvm_startup_wall_measurement": startup_measurement,
        "batches_handled": len(batches or []),
        "batch_wall_seconds": list(batches or []),
        "restarts": restarts,
        "protocol_failures": protocol_failures,
        "protocol_timeouts": protocol_timeouts,
        "shutdown_reason": shutdown_reason,
    }
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _get_project_root() -> Path:
    """Resolve project root from this file location."""
    return Path(__file__).resolve().parent.parent


def _get_predictor_dir() -> Path:
    """Return predictor directory inside the project root."""
    return _get_project_root() / "predictor"


def _get_build_dir() -> Path:
    """Return the package-local temporary Java build directory."""
    configured = os.environ.get(JAVA_BUILD_DIR_ENV, "").strip()
    build_dir = (
        Path(configured).expanduser().resolve()
        if configured
        else Path(tempfile.gettempdir()).resolve() / "demiurge_java_build"
    )
    build_dir.mkdir(parents=True, exist_ok=True)
    return build_dir


def _get_java_targets(predictor: str) -> tuple[Path, Path, str]:
    """
    Return:
    predictor_jar, batch_java, main_class
    """
    predictor_dir = _get_predictor_dir()

    if predictor == "1H":
        predictor_jar = predictor_dir / "predictorh.jar"
        batch_java = predictor_dir / "BatchProcessor1H.java"
        main_class = "predictor.BatchProcessor1H"
    elif predictor == "13C":
        predictor_jar = predictor_dir / "predictorc.jar"
        batch_java = predictor_dir / "BatchProcessor13C.java"
        main_class = "predictor.BatchProcessor13C"
    else:
        raise ValueError("predictor must be '1H' or '13C'")

    return predictor_jar, batch_java, main_class


def _ensure_java_compiled(predictor: str) -> tuple[str, str]:
    """
    Ensure Java sources are compiled globally once per project.

    Returns:
    - classpath as string
    - main_class
    """
    GREEN = '\033[38;5;46m'
    RED = '\033[38;5;196m'
    ORAN = '\033[38;5;214m'
    RESET = '\033[0m'

    project_root = _get_project_root()
    predictor_dir = _get_predictor_dir()
    build_dir = _get_build_dir()

    # These files define scientific behavior.  A missing 3D builder is not an
    # acceptable silent switch to a different 2D branch in production.
    verify_predictor_artifacts(project_root)

    cp_sep = ";" if platform.system() == "Windows" else ":"

    predictor_jar, batch_java, main_class = _get_java_targets(predictor)

    cdk_jar = predictor_dir / "cdk-2.9.jar"
    cdk_builder3d = predictor_dir / "cdk-builder3d-2.9.jar"

    for p in [predictor_jar, cdk_jar, cdk_builder3d, batch_java]:
        if not p.is_file():
            raise FileNotFoundError(f"Required file not found: {p}")

    class_rel = Path(*main_class.split(".")).with_suffix(".class")
    class_file = build_dir / class_rel

    source_sha256 = hashlib.sha256(batch_java.read_bytes()).hexdigest()
    signature_file = class_file.with_suffix(".source.sha256")
    try:
        compiled_source_sha256 = signature_file.read_text(encoding="ascii").strip()
    except FileNotFoundError:
        compiled_source_sha256 = ""

    # A content signature prevents an old hardcoded .class from surviving when
    # source timestamps are preserved while copying or unpacking the package.
    needs_compile = not class_file.exists() or compiled_source_sha256 != source_sha256

    cp_compile = [str(predictor_jar), str(cdk_jar), str(build_dir)]
    cp_compile.insert(2, str(cdk_builder3d))

    cp_run = [str(predictor_jar), str(cdk_jar), str(build_dir)]
    cp_run.insert(2, str(cdk_builder3d))

    env = os.environ.copy()
    extra_tool_opts = "-Dfile.encoding=UTF-8 -Dsun.jnu.encoding=UTF-8"
    env["JAVA_TOOL_OPTIONS"] = (env.get("JAVA_TOOL_OPTIONS", "") + " " + extra_tool_opts).strip()

    if needs_compile:
        compile_command = [
            "javac",
            "-encoding",
            "UTF-8",
            "-classpath",
            cp_sep.join(cp_compile),
            "-d",
            str(build_dir),
            "-Xlint:-options",
            "-Xlint:deprecation",
            "-proc:none",
            str(batch_java),
        ]

        try:
            subprocess.run(compile_command, check=True, env=env)
            signature_file.parent.mkdir(parents=True, exist_ok=True)
            signature_file.write_text(source_sha256 + "\n", encoding="ascii")
            print(f"\nSuccessfully compiled {main_class} into {GREEN}{build_dir}{RESET}.")
        except subprocess.CalledProcessError as e:
            print(f"{RED}Failed to compile {batch_java}: {e}{RESET}")
            raise
    else:
        print(f"\nUsing existing compiled Java class for {main_class} from {GREEN}{build_dir}{RESET}.")

    return cp_sep.join(cp_run), main_class


def clear_java_invocation_history() -> None:
    _LAST_JAVA_COMMANDS.clear()


def last_java_commands() -> dict[str, tuple[str, ...]]:
    return dict(_LAST_JAVA_COMMANDS)


def _encode_protocol_value(value: object) -> str:
    return base64.urlsafe_b64encode(str(value).encode("utf-8")).decode("ascii")


class _PersistentJavaProcessor:
    """One nucleus-specific JVM owned by the current Python worker process."""

    def __init__(
        self,
        *,
        predictor: str,
        predictor_mode: str,
        java_threads: int,
        java_heap: str,
        command: list[str],
        environment: dict[str, str],
        diagnostic_paths: dict[str, Path] | None,
        timeout_seconds: float,
    ) -> None:
        self.predictor = predictor
        self.predictor_mode = predictor_mode
        self.java_threads = java_threads
        self.java_heap = java_heap
        self.command = command
        self.environment = environment
        self.diagnostic_paths = diagnostic_paths
        self.timeout_seconds = timeout_seconds
        self.process: subprocess.Popen[str] | None = None
        self.reader_thread: threading.Thread | None = None
        self.responses: queue.Queue[str | None] = queue.Queue()
        self.started_at: float | None = None
        self.startup_seconds = 0.0
        self.startup_count = 0
        self.process_lifetime_seconds = 0.0
        self.batch_sequence = 0
        self.batches: list[dict[str, object]] = []
        self.restarts = 0
        self.protocol_failures = 0
        self.protocol_timeouts = 0
        self.shutdown_reason: str | None = None

    def _read_stdout(self) -> None:
        process = self.process
        if process is None or process.stdout is None:
            self.responses.put(None)
            return
        try:
            for line in process.stdout:
                self.responses.put(line.rstrip("\r\n"))
        finally:
            self.responses.put(None)

    def _receive(
        self,
        expected_batch: str | None = None,
        *,
        invalidate_on_error: bool = True,
    ) -> list[str]:
        try:
            line = self.responses.get(timeout=self.timeout_seconds)
        except queue.Empty as error:
            self.protocol_timeouts += 1
            if invalidate_on_error:
                self._invalidate("protocol-timeout")
            raise TimeoutError(
                f"Java {self.predictor} protocol timed out after "
                f"{self.timeout_seconds:g}s"
            ) from error
        if line is None:
            returncode = self.process.poll() if self.process is not None else None
            self.protocol_failures += 1
            if invalidate_on_error:
                self._invalidate("unexpected-eof")
            raise RuntimeError(
                f"Java {self.predictor} protocol ended unexpectedly "
                f"(returncode={returncode})"
            )
        fields = line.split("\t")
        if len(fields) < 3 or fields[0] != _PROTOCOL_PREFIX:
            self.protocol_failures += 1
            if invalidate_on_error:
                self._invalidate("malformed-response")
            raise RuntimeError(f"Malformed Java {self.predictor} protocol response: {line!r}")
        if expected_batch is not None and (
            len(fields) < 3 or fields[2] != expected_batch
        ):
            self.protocol_failures += 1
            if invalidate_on_error:
                self._invalidate("mismatched-batch-id")
            raise RuntimeError(
                f"Java {self.predictor} response batch id mismatch: {line!r}"
            )
        return fields

    def start(self) -> None:
        if self.process is not None and self.process.poll() is None:
            return
        if self.startup_count:
            self.restarts += 1
        self.responses = queue.Queue()
        started = time.perf_counter()
        self.process = subprocess.Popen(
            self.command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=None,
            text=True,
            encoding="utf-8",
            bufsize=1,
            env=self.environment,
        )
        self.started_at = time.perf_counter()
        self.startup_count += 1
        self.reader_thread = threading.Thread(
            target=self._read_stdout,
            name=f"spectraprints-java-{self.predictor}-protocol",
            daemon=True,
        )
        self.reader_thread.start()
        fields = self._receive()
        if fields != [_PROTOCOL_PREFIX, "READY", _PROTOCOL_VERSION, self.predictor]:
            self.protocol_failures += 1
            self._invalidate("invalid-ready")
            raise RuntimeError(f"Invalid Java {self.predictor} READY response: {fields!r}")
        self.startup_seconds += time.perf_counter() - started
        self._write_diagnostics()

    def run_batch(self, mol_directory: Path, output_directory: Path) -> str:
        self.start()
        process = self.process
        if process is None or process.stdin is None:
            raise RuntimeError(f"Java {self.predictor} process has no stdin")
        self.batch_sequence += 1
        batch_id = f"b{self.batch_sequence:08d}"
        request = "\t".join((
            _PROTOCOL_PREFIX,
            "BATCH",
            batch_id,
            _encode_protocol_value(mol_directory.resolve()),
            _encode_protocol_value(output_directory.resolve()),
        ))
        started = time.perf_counter()
        try:
            process.stdin.write(request + "\n")
            process.stdin.flush()
        except (BrokenPipeError, OSError) as error:
            self.protocol_failures += 1
            self._invalidate("broken-pipe")
            raise RuntimeError(f"Broken pipe to Java {self.predictor} process") from error
        fields = self._receive(batch_id)
        wall = time.perf_counter() - started
        if fields[1] == "ERROR":
            self._write_diagnostics()
            message = "unknown Java batch error"
            if len(fields) >= 4:
                try:
                    message = base64.urlsafe_b64decode(fields[3]).decode("utf-8")
                except Exception:
                    message = fields[3]
            raise RuntimeError(f"Java {self.predictor} batch {batch_id} failed: {message}")
        if fields[1] != "DONE" or len(fields) != 8:
            self.protocol_failures += 1
            self._invalidate("malformed-completion")
            raise RuntimeError(f"Malformed Java {self.predictor} completion: {fields!r}")
        try:
            record = {
                "batch_id": batch_id,
                "wall_seconds": wall,
                "input_files": int(fields[3]),
                "processed_files": int(fields[4]),
                "3d_native": int(fields[5]),
                "3d_rebuilt": int(fields[6]),
                "2d": int(fields[7]),
            }
        except ValueError as error:
            self.protocol_failures += 1
            self._invalidate("malformed-completion")
            raise RuntimeError(
                f"Malformed Java {self.predictor} completion counters: {fields!r}"
            ) from error
        self.batches.append(record)
        self._write_diagnostics()
        return str(output_directory)

    def _invalidate(self, reason: str) -> None:
        self.shutdown(reason, graceful=False)

    def shutdown(self, reason: str = "normal", graceful: bool = True) -> None:
        process = self.process
        if process is None:
            if self.shutdown_reason is None:
                self.shutdown_reason = reason
            self._write_diagnostics()
            return
        if graceful and process.poll() is None and process.stdin is not None:
            try:
                process.stdin.write(f"{_PROTOCOL_PREFIX}\tSHUTDOWN\t{_PROTOCOL_VERSION}\n")
                process.stdin.flush()
                fields = self._receive(invalidate_on_error=False)
                if fields[:3] != [_PROTOCOL_PREFIX, "BYE", _PROTOCOL_VERSION]:
                    self.protocol_failures += 1
                    raise RuntimeError(f"Invalid Java shutdown response: {fields!r}")
            except Exception:
                graceful = False
        if process.poll() is None:
            if graceful:
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    graceful = False
            if not graceful and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=5)
        for stream_name in ("stdin", "stdout"):
            stream = getattr(process, stream_name, None)
            close = getattr(stream, "close", None)
            if callable(close):
                close()
        if self.reader_thread is not None and self.reader_thread is not threading.current_thread():
            self.reader_thread.join(timeout=2)
        self.reader_thread = None
        if self.started_at is not None:
            self.process_lifetime_seconds += time.perf_counter() - self.started_at
        self.shutdown_reason = reason
        self.process = None
        self.started_at = None
        self._write_diagnostics()

    def _write_diagnostics(self) -> None:
        if self.diagnostic_paths is None:
            return
        elapsed = self.process_lifetime_seconds
        if self.started_at is not None:
            elapsed += time.perf_counter() - self.started_at
        try:
            _write_launcher_diagnostics(
                self.diagnostic_paths["launcher"],
                predictor=self.predictor,
                predictor_mode=self.predictor_mode,
                java_threads=self.java_threads,
                java_heap=self.java_heap,
                command=self.command,
                elapsed_seconds=elapsed,
                returncode=(self.process.poll() if self.process is not None else None),
                java_metrics_path=self.diagnostic_paths["java"],
                lifecycle=JAVA_LIFECYCLE_PERSISTENT,
                startup_seconds=self.startup_seconds,
                batches=self.batches,
                restarts=self.restarts,
                protocol_failures=self.protocol_failures,
                protocol_timeouts=self.protocol_timeouts,
                shutdown_reason=self.shutdown_reason,
                jvm_startup_count=self.startup_count,
                startup_measurement="popen-to-ready-handshake",
            )
        except OSError as error:
            print(f"Warning: could not write Java launcher diagnostics: {error}")


def shutdown_persistent_java_processors(reason: str = "normal") -> None:
    with _PERSISTENT_LOCK:
        processors = list(_PERSISTENT_PROCESSORS.values())
        _PERSISTENT_PROCESSORS.clear()
    for processor in processors:
        processor.shutdown(reason, graceful=not reason.startswith("python-signal-"))


atexit.register(shutdown_persistent_java_processors, "python-exit")


def run_java_batch_processor(
    mol_directory,
    predictor,
    java_threads=8,
    java_heap=DEFAULT_JAVA_HEAP,
    output_directory=None,
):
    """
    Runs the Java BatchProcessor for NMR spectrum prediction
    on the specified directory containing .mol files.

    Parameters:
    - mol_directory (str): Path to the input directory containing .mol files.
    - predictor (str): Type of NMR predictor ('1H' or '13C') to use.
    - java_threads (int): Positive Java worker count; defaults to 8 for
                          compatibility with the historical launcher.
    - java_heap (str): Positive JVM maximum heap using M/G syntax; defaults to
                       the validated production allocation (4G).

    Returns:
    - csv_output_folder (str): Path to the directory where the predicted CSV
                               files are stored.
    """

    GREEN = '\033[38;5;46m'
    RED = '\033[38;5;196m'
    ORAN = '\033[38;5;214m'
    RESET = '\033[0m'

    if isinstance(java_threads, bool) or not isinstance(java_threads, int) or java_threads <= 0:
        raise ValueError("java_threads must be a positive integer")
    java_heap = normalize_java_heap(java_heap)
    predictor_mode = _predictor_mode_from_environment()
    lifecycle = java_lifecycle_from_environment()

    csv_output_folder = (
        os.fspath(Path(output_directory).expanduser().resolve())
        if output_directory is not None
        else os.path.join(os.getcwd(), f"predicted_spectra_{predictor}")
    )
    if not os.path.exists(csv_output_folder):
        os.makedirs(csv_output_folder)
        print(f"\nCreated directory: {ORAN}{csv_output_folder}{RESET}")

    try:
        classpath, main_class = _ensure_java_compiled(predictor)
    except Exception:
        return None

    env = os.environ.copy()
    extra_tool_opts = "-Dfile.encoding=UTF-8 -Dsun.jnu.encoding=UTF-8"
    env["JAVA_TOOL_OPTIONS"] = (env.get("JAVA_TOOL_OPTIONS", "") + " " + extra_tool_opts).strip()

    print(f"\nSpectra prediction in progress (Java lifecycle={lifecycle})...\n")

    if lifecycle == JAVA_LIFECYCLE_PERSISTENT:
        key = (predictor, java_threads, java_heap, predictor_mode)
        with _PERSISTENT_LOCK:
            processor = _PERSISTENT_PROCESSORS.get(key)
            if processor is None:
                diagnostic_paths = _diagnostic_paths(predictor)
                java_command = [
                    "java",
                    f"-Xmx{java_heap}",
                    "-Dfile.encoding=UTF-8",
                    "-Dsun.jnu.encoding=UTF-8",
                ]
                if diagnostic_paths is not None:
                    java_command.append(
                        "-Xlog:gc*,safepoint:"
                        f"file={diagnostic_paths['gc']}:time,uptime,level,tags"
                    )
                    if "jfr" in diagnostic_paths:
                        java_command.append(
                            "-XX:StartFlightRecording="
                            f"filename={diagnostic_paths['jfr']},settings=profile,dumponexit=true"
                        )
                java_command.extend([
                    "-classpath",
                    classpath,
                    main_class,
                    "--persistent-service",
                    "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)",
                    "use3d",
                    str(java_threads),
                    str(diagnostic_paths["java"]) if diagnostic_paths is not None else "-",
                    predictor_mode,
                ])
                run_command = java_command
                resource_timer = Path("/usr/bin/time")
                if diagnostic_paths is not None and resource_timer.is_file():
                    run_command = [
                        str(resource_timer), "-v", "-o",
                        str(diagnostic_paths["resources"]), *java_command,
                    ]
                processor = _PersistentJavaProcessor(
                    predictor=predictor,
                    predictor_mode=predictor_mode,
                    java_threads=java_threads,
                    java_heap=java_heap,
                    command=run_command,
                    environment=env,
                    diagnostic_paths=diagnostic_paths,
                    timeout_seconds=_protocol_timeout_from_environment(),
                )
                _PERSISTENT_PROCESSORS[key] = processor
            _LAST_JAVA_COMMANDS[predictor] = tuple(processor.command)
        try:
            return processor.run_batch(Path(mol_directory), Path(csv_output_folder))
        except (OSError, RuntimeError, TimeoutError, ValueError) as error:
            print(f"{RED}Failed persistent Java {predictor} batch: {error}{RESET}")
            return None

    diagnostic_paths = _diagnostic_paths(predictor)

    java_command = [
        "java",
        f"-Xmx{java_heap}",
        "-Dfile.encoding=UTF-8",
        "-Dsun.jnu.encoding=UTF-8",
    ]

    if diagnostic_paths is not None:
        java_command.append(
            "-Xlog:gc*,safepoint:"
            f"file={diagnostic_paths['gc']}:time,uptime,level,tags"
        )
        if "jfr" in diagnostic_paths:
            java_command.append(
                "-XX:StartFlightRecording="
                f"filename={diagnostic_paths['jfr']},settings=profile,dumponexit=true"
            )

    java_command.extend([
        "-classpath",
        classpath,
        main_class,
        str(mol_directory),
        str(csv_output_folder),
        "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)",
        "use3d",
        str(java_threads),
    ])

    if diagnostic_paths is not None:
        java_command.append(str(diagnostic_paths["java"]))
    elif predictor_mode != PREDICTOR_MODE_THREAD_LOCAL:
        java_command.append("-")
    if predictor_mode != PREDICTOR_MODE_THREAD_LOCAL:
        java_command.append(predictor_mode)

    run_command = java_command
    resource_timer = Path("/usr/bin/time")
    if diagnostic_paths is not None and resource_timer.is_file():
        run_command = [
            str(resource_timer),
            "-v",
            "-o",
            str(diagnostic_paths["resources"]),
            *java_command,
        ]
    _LAST_JAVA_COMMANDS[predictor] = tuple(run_command)

    started = time.perf_counter()
    returncode = None
    try:
        completed = subprocess.run(run_command, check=True, env=env)
        returncode = completed.returncode
    except subprocess.CalledProcessError as e:
        returncode = e.returncode
        print(f"{RED}Failed to run {main_class}: {e}{RESET}")
        return None
    finally:
        if diagnostic_paths is not None:
            try:
                elapsed_seconds = time.perf_counter() - started
                startup_seconds = None
                try:
                    java_document = json.loads(
                        diagnostic_paths["java"].read_text(encoding="utf-8")
                    )
                    startup_seconds = max(
                        0.0,
                        elapsed_seconds - float(java_document["main_wall_seconds"]),
                    )
                except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
                    pass
                _write_launcher_diagnostics(
                    diagnostic_paths["launcher"],
                    predictor=predictor,
                    predictor_mode=predictor_mode,
                    java_threads=java_threads,
                    java_heap=java_heap,
                    command=run_command,
                    elapsed_seconds=elapsed_seconds,
                    returncode=returncode,
                    java_metrics_path=diagnostic_paths["java"],
                    lifecycle=JAVA_LIFECYCLE_PER_BATCH,
                    startup_seconds=startup_seconds,
                    startup_measurement="subprocess-wall-minus-java-main-wall",
                    shutdown_reason="process-exit",
                )
            except OSError as error:
                print(
                    f"{ORAN}Warning: could not write Java launcher diagnostics: "
                    f"{error}{RESET}"
                )

    return csv_output_folder
