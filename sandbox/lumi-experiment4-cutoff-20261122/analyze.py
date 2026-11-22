#!/usr/bin/env python3
"""Local immutable-ID comparison of matched old/cutoff Lumi Experiment 4.

Do not reconnect missing particles or extrapolate the aborted old run.
The CPU comparison has the same deck but not necessarily the same sampled
primary realization across execution spaces; it is not a bitwise GPU gate.
"""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DECK = "gamma_gamma_pairs-20mm-dynamic-xy-staged.in"
IP = .17
C = 299792458.
OUTPUT = HERE / "comparison"


def load_audit():
    path = ROOT / "sandbox/experiment-4/scripts/audit_loss_trajectories.py"
    spec = importlib.util.spec_from_file_location("audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def verify(run, binary_hash):
    for line in (run / "inputs.sha256").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if hashlib.sha256((run / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input mismatch {run / name}")
    if binary_hash not in (run / "execution-fingerprint.txt").read_text():
        raise ValueError("Wrong executed binary")
    log = (run / "opalx.log").read_text()
    completed = (run / "completed").exists()
    if not completed and "CIC stencil lies outside" not in log:
        raise ValueError(f"Unexpected incomplete run: {run}")
    manifest = dict(status="completed" if completed else "aborted",
                    git_revision="e24dc92185198c5f65718e057400e6c6398dce56",
                    executable_sha256=binary_hash)
    (run / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    time = re.search(r"^real ([0-9.]+)$", log, re.M)
    return dict(completed=completed, wall_seconds=float(time[1]) if time else None,
                binary_sha256=binary_hash)


def compare_prefix(left, right):
    """Compare common saved steps only; report population differences honestly."""
    rows = []
    for c in (1, 2):
        with h5py.File(left / f"{DECK[:-3]}_c{c}.h5") as a, h5py.File(right / f"{DECK[:-3]}_c{c}.h5") as b:
            for key in sorted(set(a) & set(b), key=lambda x: int(x.split("#")[-1])):
                ga, gb = a[key], b[key]
                if not np.array_equal(ga.attrs["TIME"], gb.attrs["TIME"]) or not np.array_equal(ga.attrs["GlobalTrackStep"], gb.attrs["GlobalTrackStep"]):
                    raise ValueError("Snapshot clock mismatch")
                ia, ib = ga["id"][:], gb["id"][:]
                ids, aa, bb = np.intersect1d(ia, ib, return_indices=True)
                if len(np.unique(ia)) != len(ia) or len(np.unique(ib)) != len(ib):
                    raise ValueError("Duplicate IDs")
                ra = np.column_stack([ga[k][:][aa] for k in ("x", "y", "z")]) + np.asarray(ga.attrs["RefPartR"]).reshape(3)
                rb = np.column_stack([gb[k][:][bb] for k in ("x", "y", "z")]) + np.asarray(gb.attrs["RefPartR"]).reshape(3)
                pa = np.column_stack([ga[k][:][aa] for k in ("px", "py", "pz")])
                pb = np.column_stack([gb[k][:][bb] for k in ("px", "py", "pz")])
                if not all(np.isfinite(v).all() for v in (ra, rb, pa, pb)):
                    raise ValueError("Nonfinite trajectory")
                rows.append(dict(container=c, step=int(np.asarray(ga.attrs["GlobalTrackStep"]).item()),
                                 common_ids=len(ids), only_left=len(ia)-len(ids), only_right=len(ib)-len(ids),
                                 max_delta_r_m=float(np.linalg.norm(rb-ra, axis=1).max(initial=0)),
                                 max_delta_p_mc=float(np.linalg.norm(pb-pa, axis=1).max(initial=0))))
    return pd.DataFrame(rows)


def trajectory_sample(run, c, ids=None):
    rows = []
    with h5py.File(run / f"{DECK[:-3]}_c{c}.h5") as h5:
        groups = sorted(h5.values(), key=lambda g: int(g.name.split("#")[-1]))
        for j, g in enumerate(groups):
            t = float(np.asarray(g.attrs["TIME"]).item())
            near = abs(t-IP/C) <= 20e-12
            if ids is None and j % 20 and not (near and j % 10 == 0):
                continue
            ii = g["id"][:]
            take = np.ones(len(ii), dtype=bool) if ids is None else np.isin(ii, ids)
            ref = np.asarray(g.attrs["RefPartR"]).reshape(3)
            rows.extend(zip(ii[take], np.full(take.sum(), t), g["x"][:][take]+ref[0],
                            g["z"][:][take]+ref[2], g["px"][:][take]))
    return pd.DataFrame(rows, columns=["id", "time_s", "x_m", "s_m", "px_mc"])


def main():
    OUTPUT.mkdir(exist_ok=True)
    audit, audit_path = load_audit()
    hashes = dict(old="31e97e63d5208ef210c4188b717b48ccd739fba37874f38735e08754372d1307",
                  cutoff="f44bf28796a8134cc96d8a3afbb5e8a62a90adf1be72c152dd28123266da73c3")
    summaries = {}
    for case, value in hashes.items():
        run = HERE / case
        execution = verify(run, value)
        subprocess.run([sys.executable, str(audit_path), str(run), "--ranks", "1"], check=True)
        report_path = run / "loss-audit/summary.json"
        report = json.loads(report_path.read_text())
        report["provenance"]["launch"] = "srun -n1 --gpus-per-task=1 SELECT_GPU OPALX DECK --info 4; see run-local submit.sh"
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        summaries[case] = dict(execution=execution, audit=report)
    if (HERE / "old/inputs.sha256").read_bytes() != (HERE / "cutoff/inputs.sha256").read_bytes():
        raise ValueError("Different old/new inputs")
    for name, left, right in (
        ("old_gpu_vs_cutoff_gpu", HERE / "old", HERE / "cutoff"),
        ("previous_cpu_cutoff_vs_gpu", ROOT / "sandbox/experiment-4/results/cutoff-20261122-n1", HERE / "cutoff"),
        ("previous_cpu_old_vs_gpu_old", ROOT / "sandbox/experiment-4/results/issue612-20261121-n1", HERE / "old")):
        frame = compare_prefix(left, right)
        frame.to_csv(OUTPUT / f"{name}.csv", index=False)
        summaries[name] = dict(common_snapshots=len(frame), max_delta_r_m=float(frame.max_delta_r_m.max()),
                               max_delta_p_mc=float(frame.max_delta_p_mc.max()),
                               unmatched_population_snapshots=int(((frame.only_left > 0) | (frame.only_right > 0)).sum()),
                               scope="Common saved steps and immutable IDs only; no extrapolation")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for col, case in enumerate(("old", "cutoff")):
        for c, color in ((1, "#b51f2e"), (2, "#1f5fb5")):
            tracks = trajectory_sample(HERE / case, c)
            for _, track in tracks.groupby("id"):
                axes[0, col].plot((track.s_m-IP)*1e3, track.x_m*1e3, color=color, lw=.4, alpha=.12)
                if c == 2:
                    near = abs(track.time_s-IP/C) <= 20e-12
                    axes[1, col].plot((track.time_s[near]-IP/C)*1e12, track.x_m[near]*1e6,
                                      color=color, lw=.4, alpha=.16)
        suffix = " (aborted)" if not summaries[case]["execution"]["completed"] else " (completed)"
        axes[0, col].set(title=case + suffix, xlabel=r"$s-s_{IP}$ [mm]", ylabel="x [mm]")
        axes[1, col].set(title="All positrons near nominal overlap", xlabel="Time from nominal overlap [ps]", ylabel=r"x [$\mu$m]")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Lumi Experiment 4: evolving source, 1297 witnesses/species; red electrons, blue positrons")
    fig.savefig(OUTPUT / "all_trajectories.png", dpi=200)
    plt.close(fig)
    particles = pd.read_csv(HERE / "cutoff/loss-audit/positron_identity_summary.csv")
    selected = particles.sort_values(["px_reversals", "x_crossings", "id"], ascending=[False, False, True]).head(6).id.tolist()
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), constrained_layout=True)
    for case, style, color in (("old", "--", "#d8842f"), ("cutoff", "-", "#1f5fb5")):
        tracks = trajectory_sample(HERE / case, 2, selected)
        tracks.to_csv(OUTPUT / f"{case}_selected_positrons.csv", index=False)
        for ax, i in zip(axes.flat, selected):
            track = tracks[(tracks.id == i) & (abs(tracks.time_s-IP/C) <= 20e-12)]
            ax.plot((track.time_s-IP/C)*1e12, track.x_m*1e6, ls=style, color=color, lw=1.3, label=case)
            n = int(particles.loc[particles.id == i, "px_reversals"].iloc[0])
            ax.set(title=f"ID {i}: {n} px reversals (cutoff)", xlabel="Time from nominal overlap [ps]", ylabel=r"x [$\mu$m]")
            ax.grid(alpha=.2)
    axes.flat[0].legend(frameon=False)
    fig.suptitle("Fixed-ID positron histories: six highest reversal/crossing counts; no trajectory reconnection")
    fig.savefig(OUTPUT / "positron_oscillations.png", dpi=220)
    plt.close(fig)
    (OUTPUT / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
