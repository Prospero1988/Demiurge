import os
import subprocess
import platform

def run_java_batch_processor(mol_directory, predictor):
    """
    Compiles and runs the Java BatchProcessor for NMR spectrum prediction
    on the specified directory containing .mol files.

    Parameters:
    - mol_directory (str): Path to the input directory containing .mol files.
    - predictor (str): Type of NMR predictor ('1H' or '13C') to use.

    Returns:
    - csv_output_folder (str): Path to the directory where the predicted CSV
                               files are stored.
    """

    # ANSI colors
    GREEN = '\033[38;5;46m'
    RED   = '\033[38;5;196m'
    ORAN  = '\033[38;5;214m'
    RESET = '\033[0m'

    csv_output_folder = os.path.join(os.getcwd(), f"predicted_spectra_{predictor}")
    if not os.path.exists(csv_output_folder):
        os.makedirs(csv_output_folder)
        print(f"\nCreated directory: {ORAN}{csv_output_folder}{RESET}")

    # Platform-specific classpath separator
    cp_sep = ";" if platform.system() == "Windows" else ":"

    # Common JARs
    predictor_dir = os.path.join(".", "predictor")
    cdk_jar       = os.path.join(predictor_dir, "cdk-2.9.jar")
    # Add the 3D builder JAR explicitly
    cdk_builder3d = os.path.join(predictor_dir, "cdk-builder3d-2.9.jar")

    if predictor == "1H":
        predictor_jar = os.path.join(predictor_dir, "predictorh.jar")
        batch_java    = os.path.join(predictor_dir, "BatchProcessor1H.java")
        main_class    = "predictor.BatchProcessor1H"
    elif predictor == "13C":
        predictor_jar = os.path.join(predictor_dir, "predictorc.jar")
        batch_java    = os.path.join(predictor_dir, "BatchProcessor13C.java")
        main_class    = "predictor.BatchProcessor13C"
    else:
        raise ValueError("predictor must be '1H' or '13C'")

    # Sanity checks (fail fast with a clear message)
    for p in [predictor_jar, cdk_jar]:
        if not os.path.isfile(p):
            raise FileNotFoundError(f"Required JAR not found: {p}")
    if not os.path.isfile(cdk_builder3d):
        # Not fatal — we just warn. Rebuild 3D will silently skip.
        print(f"{ORAN}Warning: cdk-builder3d jar not found at: {cdk_builder3d}{RESET}")
    if not os.path.isfile(batch_java):
        raise FileNotFoundError(f"Java source not found: {batch_java}")

    # Build classpaths
    # Compilation: builder3d is optional (reflection), but harmless to include if present
    cp_compile = [predictor_jar, cdk_jar, "."]
    if os.path.isfile(cdk_builder3d):
        cp_compile.insert(2, cdk_builder3d)   # add before "."

    # Runtime: include builder3d if present (needed for actual 3D rebuild)
    cp_run = [predictor_jar, cdk_jar, "."]
    if os.path.isfile(cdk_builder3d):
        cp_run.insert(2, cdk_builder3d)

    # --- UTF-8 fixes for JVM ---
    # 1) JVM opts wymuszające UTF-8 (stdout/stderr, parsowanie znaków, itp.)
    jvm_opts = "-Xmx8g -Dfile.encoding=UTF-8 -Dsun.jnu.encoding=UTF-8"
    # 2) Dodatkowo podajemy JAVA_TOOL_OPTIONS (JVM i tak to czyta), gdyby ktoś nadpisał file.encoding
    env = os.environ.copy()
    extra_tool_opts = "-Dfile.encoding=UTF-8 -Dsun.jnu.encoding=UTF-8"
    env["JAVA_TOOL_OPTIONS"] = (env.get("JAVA_TOOL_OPTIONS", "") + " " + extra_tool_opts).strip()

    # (opcjonalnie) na Windows ustaw konsolę na UTF-8, żeby inne narzędzia też drukowały poprawnie
    if platform.system() == "Windows":
        try:
            subprocess.run("chcp 65001 >NUL", shell=True, check=False)
        except Exception:
            pass

    compile_command = (
        f'javac -classpath "{cp_sep.join(cp_compile)}" '
        f'-d . -Xlint:-options -Xlint:deprecation -proc:none "{batch_java}"'
    )

    try:
        subprocess.run(compile_command, shell=True, check=True)
        print(f"\nSuccessfully compiled {main_class}.")
        print("\nSpectra prediction in progress...\n")
    except subprocess.CalledProcessError as e:
        print(f"{RED}Failed to compile {batch_java}: {e}{RESET}")
        return None

    # JVM options: feel free to tweak heap if needed
    jvm_opts = "-Xmx8g"

    # Run
    run_command = (
        f'java {jvm_opts} -classpath "{cp_sep.join(cp_run)}" '
        f'{main_class} "{mol_directory}" "{csv_output_folder}" '
        f'"Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)"'
    )

    try:
        subprocess.run(run_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"{RED}Failed to run {main_class}: {e}{RESET}")
        return None

    return csv_output_folder
