#!/usr/bin/env python3
"""Record a DBA amplitude/Richardson experiment using runtime OPTION settings.

This runner never modifies production source or the user's input template.
Use a separate output directory for each Richardson level.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

from convergence_dt import DBA, plt, pd, np, run_case

REPOSITORY = DBA.parents[2]
BUILDER = REPOSITORY / "src/Algorithms/LinearTransferMapBuilder.cpp"


def configure_map(input_text: str, epsilon: float, levels: int, integrator: str = "BORIS") -> str:
    if integrator not in ("BORIS", "RK4", "DOP853"):
        raise ValueError("unsupported map integrator")
    steps = ",".join([f"{epsilon:.16e}"] * 6)
    option = (f"OPTION, ENABLELINEARTRANSFERMAPS=TRUE, LINEARTRANSFERMAPRICHARDSON={levels}, "
              f'LINEARTRANSFERMAPSTEPS={{{steps}}}, LINEARTRANSFERMAPINTEGRATOR="{integrator}";\n\n')
    updated, count = re.subn(r"(?m)^TRACK\b", option + "TRACK", input_text)
    if count != 1:
        raise RuntimeError("expected exactly one TRACK statement")
    return updated


def summarize(output: Path) -> pd.DataFrame:
    rows = []
    for result_file in output.glob("eps-*/ranks-*/result.json"):
        result = json.loads(result_file.read_text())
        result.setdefault("richardson_levels", 0)  # Original pre-option study provenance.
        rank_directories = sorted(result_file.parent.parent.glob("ranks-*"),
                                  key=lambda path: int(path.name.split("-")[1]))
        matrix_name = "map-2-dba-dt-1e-13-matrix.txt"
        reference = np.loadtxt(rank_directories[0] / matrix_name)
        measured = np.loadtxt(result_file.parent / matrix_name)
        result["mpi_matrix_max_difference"] = float(np.max(np.abs(measured - reference)))
        rows.append({key: value for key, value in result.items()
                     if key not in ("provenance", "steps")})
    frame = pd.DataFrame(rows).sort_values(["ranks", "epsilon"], ascending=[True, False])
    if frame.richardson_levels.nunique() != 1:
        raise RuntimeError("use a separate output directory for each Richardson level")
    levels = int(frame.richardson_levels.iloc[0])
    frame.to_csv(output / "convergence-perturbations.csv", index=False, float_format="%.12e")
    figure, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), constrained_layout=True)
    for axis, column, label in zip(axes, ("abs_R16", "abs_R26"),
                                   (r"$|R_{16}|$ [m]", r"$|R_{26}|$")):
        for ranks, group in frame.groupby("ranks"):
            axis.loglog(group.epsilon, group[column], marker="o" if ranks == 1 else "x",
                        linestyle="-" if ranks == 1 else "--", linewidth=1.5,
                        label=f"{ranks} MPI rank{'s' if ranks != 1 else ''}")
        baseline = frame[frame.ranks == frame.ranks.min()].iloc[0]
        # Show the asymptotic guide near the coarsest points without compressing
        # the measured residuals to make room for an unobserved extrapolation.
        amplitudes = np.sort(frame.epsilon.unique())[-2:]
        order = 2 * (levels + 1)
        axis.loglog(amplitudes, baseline[column] * (amplitudes / baseline.epsilon) ** order,
                    ":", color="0.5", label=rf"$\propto\epsilon^{{{order}}}$ guide")
        axis.set_xlabel(r"Coordinate perturbation $\epsilon$ (m or dimensionless)")
        axis.set_ylabel(label)
        axis.grid(True, which="major", alpha=0.25)
        axis.grid(True, which="minor", alpha=0.08)
        axis.invert_xaxis()
    axes[0].legend(frameon=False, fontsize=9)
    figure.suptitle(r"DBA dispersion residuals; $\Delta t=10^{-13}$ s" +
                    f"; Richardson levels={levels}")
    figure.savefig(output / "convergence-perturbations.png", dpi=220)
    plt.close(figure)
    columns = ["epsilon", "richardson_levels", "ranks", "R16", "R26", "full_matrix_error",
               "determinant_error", "canonical_J_error", "mpi_matrix_max_difference"]
    print(frame[columns].to_string(index=False, float_format=lambda value: f"{value:.8e}"))
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=REPOSITORY / "omp-build")
    parser.add_argument("--input", type=Path, default=DBA / "map-2-dba.in")
    parser.add_argument("--output-dir", type=Path, default=DBA / "perturbations")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1],
                        help="MPI ranks per case; default: single rank")
    parser.add_argument("--mpi-args", default="--map-by slot:OVERSUBSCRIBE --bind-to none")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--epsilon", type=float, default=1.e-3)
    parser.add_argument("--richardson-levels", type=int, choices=range(5), default=0)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    if args.summarize_only:
        summarize(args.output_dir)
        return 0
    if any(ranks < 1 for ranks in args.ranks) or args.jobs < 1:
        parser.error("ranks and jobs must be positive")
    if not np.isfinite(args.epsilon) or not 0 < args.epsilon < 1:
        parser.error("--epsilon must be finite and strictly between zero and one")
    for existing in args.output_dir.glob("eps-*/ranks-*/result.json"):
        if json.loads(existing.read_text()).get("richardson_levels", 0) != args.richardson_levels:
            parser.error("use a separate --output-dir for each Richardson level")
    source = BUILDER.read_text()
    steps = [args.epsilon] * 6
    input_text = args.input.read_text()
    if re.search(r"\bPSDUMPFREQ\s*=\s*[,;]", input_text):
        parser.error("the template has an empty PSDUMPFREQ; supply a valid --input copy")
    input_text = configure_map(input_text, args.epsilon, args.richardson_levels)
    mantissa, exponent = f"{steps[0]:.8e}".split("e")
    amplitude_label = mantissa.rstrip("0").rstrip(".") + "e" + exponent
    amplitude_dir = args.output_dir / f"eps-{amplitude_label}"
    amplitude_dir.mkdir(parents=True, exist_ok=True)
    with (amplitude_dir / "build.log").open("w") as build_log:
        subprocess.run(["cmake", "--build", str(args.build_dir), "--target", "opalx_exe",
                        "-j", str(args.jobs)], stdout=build_log, stderr=subprocess.STDOUT,
                       check=True)
    if BUILDER.read_text() != source:
        raise RuntimeError("builder source changed during the build")
    executable = (args.build_dir / "src/opalx").resolve()
    provenance = {
        "git_revision": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                 cwd=REPOSITORY, text=True).strip(),
        "builder_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "input_template": input_text,
        "input_path": str(args.input.resolve()),
        "mpi_args": args.mpi_args,
    }
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OMP_PROC_BIND", "false")
    for ranks in args.ranks:
        output = amplitude_dir / f"ranks-{ranks}"
        output.mkdir(parents=True, exist_ok=True)
        print(f"Running epsilon={steps[0]:.8g}, Richardson={args.richardson_levels}, "
              f"DT=1e-13 s, ranks={ranks}", flush=True)
        with tempfile.TemporaryDirectory(prefix="opalx-map-perturbations-") as directory:
            result = run_case(executable, input_text, 1e-13, Path(directory), output,
                              ranks, tuple(shlex.split(args.mpi_args)))
        if result["status"] != "OK":
            raise RuntimeError(f"study failed; inspect {output}")
        stdout = (output / "map-2-dba-dt-1e-13.out").read_text()
        if f"Richardson levels: {args.richardson_levels};" not in stdout:
            raise RuntimeError("stdout does not confirm the requested Richardson setting")
        result.update(epsilon=steps[0], steps=steps, richardson_levels=args.richardson_levels,
                      ranks=ranks, provenance=provenance)
        (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    summarize(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
