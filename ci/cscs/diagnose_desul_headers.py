#!/usr/bin/env python3
"""Report Desul header selection using the actual Ninja/Kokkos compiler command.

This is a preprocessing-only CI probe: it does not compile objects, change include
precedence, or execute OPALX. Its caller must preserve the normal build result.
"""

import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import tempfile


OBJECT_TARGET = "src/CMakeFiles/opalx.dir/BeamlineCore/RFCavityRep.cpp.o"


def preprocess_arguments(command):
    """Keep Ninja's launchers and flags, removing object/dependency generation."""
    args = shlex.split(command)
    if "-c" not in args or any(arg in {"&&", "||", ";", "|"} for arg in args):
        raise ValueError("Expected a single compile command from Ninja")
    result = []
    index = 0
    while index < len(args):
        arg = args[index]
        if arg in {"-o", "-MF", "-MT", "-MQ"}:
            if index + 1 == len(args):
                raise ValueError("Missing argument after " + arg)
            index += 2
            continue
        if arg in {"-c", "-MD", "-MMD", "-MP", "-MG", "-M", "-MM"}:
            index += 1
            continue
        if arg.startswith(("-o", "-MF", "-MT", "-MQ")):
            index += 1
            continue
        result.append(arg)
        index += 1
    return result + ["-E"]


def selected_headers(preprocessed):
    """Read real include paths from preprocessor line markers, without API guesses."""
    paths = set()
    marker = re.compile(r'^\s*#\s*(?:line\s+)?\d+\s+"([^"]+)"')
    with preprocessed.open(errors="replace") as stream:
        for line in stream:
            match = marker.match(line)
            if not match:
                continue
            path = match.group(1)
            if "/desul/" in path or path.startswith("desul/") or path.endswith("Kokkos_Atomics_Desul_Wrapper.hpp"):
                paths.add(path)
    return sorted(paths)


def main(argv=None):
    """Preprocess the failing translation unit and print a bounded CI report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", type=Path)
    options = parser.parse_args(argv)
    build_dir = options.build_dir.resolve()
    print("Desul header-selection diagnostic (preprocessing only)", flush=True)
    for name in ("CC", "CXX", "CUDAHOSTCXX", "NVCC_WRAPPER_DEFAULT_COMPILER"):
        print("{}={}".format(name, os.environ.get(name, "<unset>")))
    for name in ("gcc", "g++", "nvcc"):
        print("{} on PATH: {}".format(name, shutil.which(name) or "<not found>"))

    try:
        # compile_commands.json omits RULE_LAUNCH_COMPILE: use Ninja to retain
        # kokkos_launch_compiler/nvcc_wrapper and the actual CUDA include search.
        commands = subprocess.run(
            ["ninja", "-t", "commands", OBJECT_TARGET],
            cwd=build_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=30,
        )
        lines = commands.stdout.strip().splitlines()
        if not lines:
            raise ValueError("Ninja returned no compiler command")
        args = preprocess_arguments(lines[-1])
        if not any(arg.endswith("/BeamlineCore/RFCavityRep.cpp") for arg in args):
            raise ValueError("Ninja command does not select RFCavityRep.cpp")
        print("Preprocessor command (includes launchers and ordered include flags):")
        print(" ".join(shlex.quote(arg) for arg in args), flush=True)
        with tempfile.TemporaryDirectory(prefix="opalx-desul-") as temporary:
            output = Path(temporary) / "preprocessed.ii"
            errors = Path(temporary) / "stderr.log"
            with output.open("w") as stdout, errors.open("w") as stderr:
                process = subprocess.Popen(
                    args, cwd=build_dir, stdout=stdout, stderr=stderr, start_new_session=True,
                )
                try:
                    returncode = process.wait(timeout=90)
                except subprocess.TimeoutExpired:
                    # Kokkos launches NVCC through shell wrappers. Stop the whole
                    # diagnostic process group, including child preprocessors.
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
                    raise
            print("Preprocessor exit status:", returncode)
            stderr_lines = errors.read_text(errors="replace").splitlines()
            if stderr_lines:
                print("Preprocessor stderr (last 40 lines):")
                print("\n".join(stderr_lines[-40:]))
            headers = selected_headers(output)
            print("Selected Desul/Kokkos wrapper headers:")
            for header in headers:
                print("  {}\n    resolved: {}".format(header, (build_dir / header).resolve()))
            if not headers:
                print("  No matching line markers; header selection is inconclusive.")
            return 0 if returncode == 0 and headers else 1
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print("Diagnostic incomplete:", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
