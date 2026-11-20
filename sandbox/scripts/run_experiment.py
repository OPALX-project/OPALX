#!/usr/bin/env python3
"""Local scientific experiment launcher; deliberately separate from CI regressions."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

SANDBOX = Path(__file__).resolve().parents[1]
ROOT = SANDBOX.parent
TITLES = {
    1: "Rigid two-Gaussian analytic/PIC field comparison",
    2: "Exact-timed twelve-witness CAIN comparison",
    3: "Full-population timed CAIN emission",
    4: "32 cm production geometry and field-window lifecycle",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def summarize_h5(run):
    """Report finite phase space and plotted RMS sizes, without claiming convergence."""
    import h5py
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = sorted(run.glob("*_c*.h5")) or [
        run / (deck.stem + ".h5") for deck in sorted(run.glob("*.in"))
        if (run / (deck.stem + ".h5")).is_file()
    ]
    if not files:
        raise ValueError("No particle-container H5 files produced")
    report = {}
    nonempty_containers = 0
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for file in files:
        records = []
        with h5py.File(file, "r") as h5:
            for key in sorted(h5, key=lambda k: int(k.split("#")[-1])):
                group = h5[key]
                if not all(k in group for k in ("x", "y", "z", "px", "py", "pz")):
                    continue
                arrays = [group[k][:] for k in ("x", "y", "z", "px", "py", "pz")]
                if any(not np.isfinite(a).all() for a in arrays):
                    raise ValueError(f"Non-finite phase space in {file.name}:{key}")
                time = float(np.asarray(group.attrs["TIME"]).reshape(-1)[0])
                step = int(np.asarray(group.attrs.get("GlobalTrackStep", [-1])).reshape(-1)[0])
                records.append((time, len(arrays[0]), float(np.std(arrays[0])) if len(arrays[0]) else 0., step))
        if not records:
            report[file.stem] = {"dumps": 0, "finite_phase_space": None,
                                 "note": "No particle snapshots; primary dumps may be suppressed by the input"}
            continue
        nonempty_containers += 1
        values = np.asarray(records)
        container = re.search(r"_c(\d+)$", file.stem)
        label = f"container {container[1]}" if container else "primary"
        if len(records) == 1:
            label += " (one stored snapshot)"
        marker = "o" if len(records) == 1 else None
        axes[0].plot(values[:, 0] * 1e12, values[:, 1], label=label, marker=marker)
        axes[1].plot(values[:, 0] * 1e12, values[:, 2] * 1e6, label=label, marker=marker)
        report[file.stem] = {"dumps": len(records), "final_count": records[-1][1],
                            "last_global_step": records[-1][3], "finite_phase_space": True}
    if not nonempty_containers:
        raise ValueError("No phase-space samples in any container")
    axes[0].set_ylabel("Macroparticle count")
    axes[1].set_ylabel("RMS x [µm]")
    for ax in axes:
        ax.set_xlabel("Time [ps]")
        ax.grid(alpha=.25)
        ax.legend(fontsize=7)
    fig.savefig(run / "particle_summary.png", dpi=200)
    fig.savefig(run / "particle_summary.pdf")
    plt.close(fig)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", type=int, choices=TITLES, required=True)
    parser.add_argument("--merlin6", action="store_true", help="Reserved; A100 submission is not implemented yet")
    parser.add_argument("--smoke", action="store_true", help="Quick execution check, not a converged physics result")
    parser.add_argument("--full", action="store_true", help="Experiment 4: track the full population with the staged production deck")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--analyze-only", type=Path, metavar="RUN_DIR")
    parser.add_argument("--output-dir", type=Path, help="New, empty directory below this experiment's results/")
    parser.add_argument("--opalx", type=Path, default=Path(os.environ.get("OPALX_EXE", ROOT / "build_openmp/src/opalx")))
    parser.add_argument("--threads", type=int, choices=range(1, 5), default=4)
    parser.add_argument("--particles", type=int)
    parser.add_argument("--mesh", type=int, nargs=3, metavar=("NX", "NY", "NZ"))
    args = parser.parse_args()
    if args.merlin6:
        parser.error("--merlin6 is reserved: A100 submission is not implemented; no job was submitted")
    if args.full and (args.experiment != 4 or args.smoke):
        parser.error("--full is only supported for experiment 4 and cannot be combined with --smoke")
    if args.experiment in (3, 4) and (args.particles or args.mesh):
        parser.error("--particles/--mesh apply only to experiments 1 and 2; other decks retain their parameters")
    if args.particles is not None and args.particles <= 0:
        parser.error("--particles must be positive")
    if args.mesh and min(args.mesh) < 2:
        parser.error("mesh dimensions must be at least two")
    if args.analyze_only and (args.output_dir or args.prepare_only):
        parser.error("--analyze-only cannot be combined with --output-dir/--prepare-only")
    experiment = SANDBOX / f"experiment-{args.experiment}"
    inputs, scripts = experiment / "inputs", experiment / "scripts"
    result_root = experiment / "results"
    result_root.mkdir(exist_ok=True)
    if args.analyze_only:
        run = args.analyze_only.resolve()
        if not run.is_relative_to(result_root.resolve()) or not run.is_dir():
            parser.error("analysis directory must exist below this experiment's results/")
        previous = json.loads((run / "run_manifest.json").read_text())
        if previous["experiment"] != args.experiment:
            parser.error("run belongs to a different experiment")
        args.smoke = previous["smoke"]
        args.full = previous["full"]
    else:
        if not args.prepare_only and not args.opalx.is_file():
            parser.error(f"OPALX executable not found: {args.opalx}")
        if args.output_dir:
            run = args.output_dir.resolve()
            if not run.is_relative_to(result_root.resolve()) or run == result_root.resolve():
                parser.error("--output-dir must be a new subdirectory of results/")
            if run.exists():
                parser.error("refusing to overwrite an existing output directory")
            run.mkdir(parents=True)
        else:
            run = Path(tempfile.mkdtemp(prefix=datetime.now().strftime("local-%Y%m%d-%H%M%S-"), dir=result_root))
    env = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = str(args.threads)
    env.update(MPLBACKEND="Agg", MPLCONFIGDIR=str(run / ".plot-cache"),
               XDG_CACHE_HOME=str(run / ".cache"), PYTHONDONTWRITEBYTECODE="1")
    os.environ.update(env)
    provenance = {"experiment": args.experiment, "title": TITLES[args.experiment], "smoke": args.smoke,
                  "full": args.full, "started_utc": datetime.now(timezone.utc).isoformat(), "commands": []}
    if args.analyze_only:
        provenance = previous
    else:
        provenance["git_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        provenance["git_status"] = subprocess.check_output(["git", "status", "--short", "--untracked-files=no"], cwd=ROOT, text=True)
        if args.opalx.is_file():
            provenance["executable_sha256"] = sha256(args.opalx)
        source_files = [p for folder in [SANDBOX / "scripts", *SANDBOX.glob("experiment-*/scripts")]
                        for p in folder.glob("*.py")]
        provenance["source_sha256"] = {str(p.relative_to(SANDBOX)): sha256(p) for p in
            [*source_files, *inputs.glob("*")] if p.is_file()}

    def execute(command, label):
        command = [str(x) for x in command]
        print(f"[{label}] {' '.join(command)}", flush=True)
        provenance["commands"].append(command)
        with (run / f"{label}.log").open("w") as log:
            subprocess.run(command, cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)

    def py(name, *options, label=None):
        execute([sys.executable, scripts / name, *options], label or Path(name).stem)

    def opalx(deck):
        execute([args.opalx.resolve(), run / deck, "--info", "4"], "opalx")
        if not (run / "timing.dat").is_file() or not (run / "timing.dat").stat().st_size:
            raise ValueError("OPALX did not write a nonempty timing.dat")

    def analyze():
        if args.experiment == 1:
            py("compare_opalx_fields.py", "--case-dir", run / "cases", "--output-dir", run / "analysis")
            import pandas as pd
            table = pd.read_csv(run / "analysis/opalx_vs_rigid_gaussian_summary.csv")
            metrics = table[["case", "relative_l2_E", "relative_l2_B", "uncovered_probe_fraction"]]
            if len(table) != 4 or not all(math.isfinite(float(value)) for value in metrics.iloc[:, 1:].to_numpy().flat):
                raise ValueError("Expected finite comparison metrics for all four source separations")
            return {"field_comparison": "analysis/opalx_vs_rigid_gaussian_summary.csv",
                    "field_metrics": metrics.to_dict(orient="records")}
        report = summarize_h5(run)
        if args.experiment == 2 and not args.smoke:
            py("compare_timed_track12.py", "--run-dir", run)
            report["trajectory_comparison"] = "results/track12_comparison.json"
            report["trajectory_metrics"] = json.loads((run / report["trajectory_comparison"]).read_text())["overall"]
        if args.experiment == 3:
            sys.path.insert(0, str(scripts))
            validation = load(scripts / "validate_opalx_output.py", "emission_validation")
            report["electron_births"] = validation.validate_witness(run / "gamma_gamma_pairs_timed_c1.h5", run / "input/cain_electrons.emittedfromfile", 1, -1)
            report["positron_births"] = validation.validate_witness(run / "gamma_gamma_pairs_timed_c2.h5", run / "input/cain_positrons.emittedfromfile", 2, 1)
            if report["gamma_gamma_pairs_timed_c0"]["final_count"] != 10000 or report["gamma_gamma_pairs_timed_c0"]["last_global_step"] != 159:
                raise ValueError("The 160-step full-population run did not retain its 10,000-particle primary")
        if args.experiment == 4 and not args.full:
            states = list(dict.fromkeys(re.findall(r"BB-state=(\w+)", (run / "opalx.log").read_text())))
            if states != ["Inactive", "Active", "Completed"]:
                raise ValueError(f"Incomplete field-window lifecycle: {states}")
            primary = report["gamma_gamma_pairs-20mm-dynamic-xy-smoke"]
            if primary["last_global_step"] != 1607 or primary["final_count"] != 32:
                raise ValueError("Lifecycle smoke did not finish 1,607 steps with 32 primary particles")
            report["field_window_states"] = states
        return report

    try:
        if not args.analyze_only:
            if args.experiment == 1:
                particles = args.particles or (1000 if args.smoke else 20000)
                mesh = args.mesh or ([8, 8, 16] if args.smoke else [32, 32, 64])
                if mesh[0] != mesh[1]:
                    raise ValueError("experiment 1 generator requires NX=NY")
                py("rigid_two_gaussian_fields.py", "--output-dir", run / "analytic")
                options = ["--output-dir", run / "cases", "--primary-macroparticles", particles,
                           "--nxy", mesh[0], "--nz", mesh[2], "--probe-nx", 21, "--probe-nz", 41,
                           "--omp-threads", args.threads, "--opalx", args.opalx.resolve()]
                if not args.prepare_only:
                    options.append("--run")
                py("make_opalx_field_cases.py", *options)
            elif args.experiment == 2:
                mesh = args.mesh or [16, 16, 32]
                py("prepare_timed_track12.py", "--output-dir", run, "--primary-macroparticles", args.particles or 2000,
                   "--nx", mesh[0], "--ny", mesh[1], "--nz", mesh[2])
                if args.smoke:
                    deck = run / "track12_timed.in"
                    text, count = re.subn(r"REAL\s+n_steps\s*=\s*1501", "REAL n_steps=3", deck.read_text())
                    if count != 1:
                        raise ValueError("Unexpected Track12 MAXSTEPS; smoke deck not generated")
                    deck.write_text(text)
                if not args.prepare_only:
                    opalx("track12_timed.in")
            elif args.experiment == 3:
                shutil.copy2(inputs / "gamma_gamma_pairs_timed.in", run)
                py("convert_cain_pairs.py", "--input", inputs / "fort98.txt", "--output-dir", run / "input", "--manifest", run / "conversion_manifest.json")
                py("validate_pair_input.py", "--input", inputs / "fort98.txt", "--output-dir", run / "input",
                   "--deck", run / "gamma_gamma_pairs_timed.in", "--report", run / "input_validation.json")
                if not args.prepare_only:
                    opalx("gamma_gamma_pairs_timed.in")
            else:
                name = "gamma_gamma_pairs-20mm-dynamic-xy-" + ("staged.in" if args.full else "smoke.in")
                text = (inputs / name).read_text().replace("../track-e-p/", "input/")
                (run / name).write_text(text)
                (run / "input").mkdir()
                for file in (inputs / "legacy").glob("*.fromfile"):
                    shutil.copy2(file, run / "input" / file.name)
                py("analyze_pair_momenta.py", "--electron-file", run / "input/gamma_gamma_electrons.fromfile",
                   "--positron-file", run / "input/gamma_gamma_positrons.fromfile", "--output-dir", run / "input-analysis", "--plot")
                if not args.prepare_only:
                    opalx(name)
        if args.prepare_only:
            provenance["status"] = "prepared, not simulated"
        else:
            provenance["analysis"] = analyze()
            provenance["status"] = "completed"
        provenance.pop("error", None)
        (run / "REPORT.md").write_text(
            f"# {TITLES[args.experiment]}\n\nStatus: {provenance['status']}.\n\n"
            "Local coarse/smoke settings are workflow checks, not convergence evidence. "
            "See run_manifest.json for exact inputs, commands and analysis.\n\n"
            + ("This lifecycle smoke tracks only the primary; --full enables both witness populations.\n" if args.experiment == 4 and not args.full else "")
            + ("\n## Analysis summary\n\n```json\n" + json.dumps(provenance["analysis"], indent=2) + "\n```\n" if "analysis" in provenance else "")
        )
    except Exception as error:
        provenance.update(status="failed", error=str(error))
        (run / "REPORT.md").write_text(f"# {TITLES[args.experiment]}\n\nStatus: failed.\n\n{error}\n\nSee the preserved logs and run_manifest.json.\n")
        print(f"FAILED: {error}\nPreserved logs: {run}", file=sys.stderr)
        return 1
    finally:
        provenance["updated_utc"] = datetime.now(timezone.utc).isoformat()
        provenance["input_sha256"] = {str(p.relative_to(run)): sha256(p) for p in run.rglob("*")
                                      if p.is_file() and p.suffix in (".in", ".fromfile", ".emittedfromfile")}
        (run / "run_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(f"Results: {run}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
