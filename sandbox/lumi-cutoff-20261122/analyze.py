#!/usr/bin/env python3
"""Validate a matched Lumi Figure 6(b) rerun and compare immutable H5 records.

Run locally with ~/.venv-h6 after fetching completed output. The old run and
its Figure 6(b) assets are never modified. No acceptance tolerance is inferred.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(variable, "1")

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integrity(run, ip):
    rows = []
    for species in (1, 2):
        previous = set()
        canonical = hashlib.sha256()
        records, outside, max_s, frames = 0, 0, 0.0, 0
        with h5py.File(run / f"track12_timed_c{species}.h5", "r") as h5:
            for key in sorted(h5, key=lambda name: int(name.split("#")[-1])):
                group = h5[key]
                ids = group["id"][:]
                current = set(map(int, ids))
                if len(current) != len(ids) or not previous <= current:
                    raise ValueError(f"{run}/{key}: duplicate or disappearing witness IDs")
                if any(not np.isfinite(group[k][:]).all() for k in ("x", "y", "z", "px", "py", "pz")):
                    raise ValueError(f"{run}/{key}: nonfinite phase space")
                if np.any(np.asarray(group.attrs["TaitBryantAngles"]) != 0):
                    raise ValueError("Nonzero reference rotation needs an explicit lab transformation")
                # Canonical order by immutable ID. Hash actual values/metadata
                # bytewise, rather than mistaking max_abs==0 for byte equality
                # (e.g. +0.0 and -0.0 have different bit patterns).
                order = np.argsort(ids)
                canonical.update(ids[order].tobytes())
                for component in ("x", "y", "z", "px", "py", "pz", "q", "sp"):
                    canonical.update(group[component][:][order].tobytes())
                for attribute in ("TIME", "GlobalTrackStep", "RefPartR", "TaitBryantAngles"):
                    canonical.update(np.asarray(group.attrs[attribute]).tobytes())
                relative_s = group["z"][:] + np.asarray(group.attrs["RefPartR"]).reshape(3)[2] - ip
                outside += int(np.count_nonzero(np.abs(relative_s) > .01))
                max_s = max(max_s, float(np.max(np.abs(relative_s), initial=0)))
                records += len(ids)
                frames += 1
                previous = current
                final_step = int(np.asarray(group.attrs["GlobalTrackStep"]).item())
        if final_step != 1500 or len(previous) != 6 or records != 6504:
            raise ValueError(f"Incomplete Track12 container {species}: {final_step}, {len(previous)}, {records}")
        rows.append(dict(container=species, snapshots=frames, records=records,
                         final_count=len(previous), final_global_step=final_step,
                         finite=True, ids_persist=True, outside_window_records=outside,
                         max_distance_from_ip_m=max_s, canonical_h5_sha256=canonical.hexdigest()))
    return rows


def paired_norms(left, right):
    count, max_r, max_p, sum_r, sum_p = 0, 0., 0., 0., 0.
    for species in (1, 2):
        with h5py.File(left / f"track12_timed_c{species}.h5") as a, h5py.File(right / f"track12_timed_c{species}.h5") as b:
            for key in a:
                ga, gb = a[key], b[key]
                ia, ib = np.argsort(ga["id"][:]), np.argsort(gb["id"][:])
                ra = np.column_stack([ga[k][:][ia] for k in ("x", "y", "z")]) + np.asarray(ga.attrs["RefPartR"]).reshape(3)
                rb = np.column_stack([gb[k][:][ib] for k in ("x", "y", "z")]) + np.asarray(gb.attrs["RefPartR"]).reshape(3)
                pa = np.column_stack([ga[k][:][ia] for k in ("px", "py", "pz")])
                pb = np.column_stack([gb[k][:][ib] for k in ("px", "py", "pz")])
                dr2, dp2 = np.sum((rb - ra)**2, axis=1), np.sum((pb - pa)**2, axis=1)
                count += len(ia)
                sum_r += float(dr2.sum())
                sum_p += float(dp2.sum())
                max_r = max(max_r, float(dr2.max(initial=0)))
                max_p = max(max_p, float(dp2.max(initial=0)))
    return dict(records=count, max_delta_lab_r_m=np.sqrt(max_r), max_delta_p_mc=np.sqrt(max_p),
                rms_delta_lab_r_m=np.sqrt(sum_r/count), rms_delta_p_mc=np.sqrt(sum_p/count))


def make_plot(merged, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(11, 5.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    tracks = fig.add_subplot(grid[:, 0])
    dr = fig.add_subplot(grid[0, 1])
    dp = fig.add_subplot(grid[1, 1], sharex=dr)
    for (species, pair), group in merged.groupby(["species", "pair"]):
        group = group.sort_values("reference_step")
        color = "#b51f2e" if species == "electron" else "#1f5fb5"
        tracks.plot(group.s_opalx_m_old * 1e6, group.x_opalx_m_old * 1e6,
                    color=color, lw=2.0, alpha=.35)
        tracks.plot(group.s_opalx_m_new * 1e6, group.x_opalx_m_new * 1e6,
                    color=color, lw=.9, ls="--")
        label = species if pair == 1 else None
        dr.plot(group.ct_m_new * 1e3, group.delta_r_m * 1e9, color=color, lw=.8, label=label)
        dp.plot(group.ct_m_new * 1e3, group.delta_p_mc, color=color, lw=.8)
    tracks.set(xlabel=r"$s-s_{IP}$ [$\mu$m]", ylabel=r"$x$ [$\mu$m]",
               title="Previous Figure 6(b) and cutoff rerun")
    tracks.legend(handles=[Line2D([], [], color="k", lw=2, alpha=.35, label="previous"),
                           Line2D([], [], color="k", ls="--", label="cutoff"),
                           Line2D([], [], color="#b51f2e", label="electron"),
                           Line2D([], [], color="#1f5fb5", label="positron")],
                  fontsize=8, frameon=False)
    dr.set(ylabel=r"$|\Delta\mathbf{r}|$ [nm]", title="Cutoff minus previous; matched particles and times")
    dp.set(xlabel=r"CAIN $ct$ [mm]", ylabel=r"$|\Delta\mathbf{p}|$ [$m_ec$]")
    dr.legend(fontsize=8, frameon=False)
    for ax in (tracks, dr, dp):
        ax.grid(alpha=.2)
    fig.suptitle("Lumi Track-12: one MI250, 400,000 primaries, 256×256×128, 1,501 steps")
    fig.savefig(output / "previous_vs_cutoff.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=HERE / "track12")
    parser.add_argument("--baseline", type=Path, default=ROOT / "sandbox/lumi-validation-20261120/track12")
    parser.add_argument("--control", type=Path, help="Optional repeat with the saved old executable")
    args = parser.parse_args()
    run, baseline = args.run.resolve(), args.baseline.resolve()
    if run == baseline:
        raise ValueError("Refusing to modify the baseline directory")
    if not (run / "completed").is_file() or not (run / "timing.dat").stat().st_size:
        raise ValueError("The remote case did not complete")
    preparation = json.loads((baseline / "results/preparation_manifest.json").read_text())
    for name, expected in preparation["input_sha256"].items():
        if digest(baseline / name) != expected or digest(run / name) != expected:
            raise ValueError(f"Input changed: {name}")
    if preparation["mesh"] != [256, 256, 128] or preparation["primary"]["particle_count"] != 400000:
        raise ValueError("Wrong Figure 6(b) baseline")
    model = load(ROOT / "sandbox/experiment-2/scripts/compare_timed_track12.py", "track12_comparison")
    ip = model.interaction_point(run / "track12_timed.in")
    checks = integrity(run, ip)
    old_checks = integrity(baseline, ip)
    comparer = load(ROOT / "sandbox/scripts/compare_smoke_phase_space.py", "raw_h5_comparison")
    differences = comparer.compare(baseline, run)
    output = run / "cutoff-comparison"
    output.mkdir(exist_ok=True)
    differences.to_csv(output / "h5_differences.csv", index=False)
    stats = differences.groupby(["file", "component"]).agg(
        max_abs=("max_abs", "max"), squared_error=("squared_error", "sum"), count=("count", "sum"))
    stats["rms"] = np.sqrt(stats.squared_error / stats["count"])
    stats.to_csv(output / "h5_summary.csv")
    subprocess.run([sys.executable, str(ROOT / "sandbox/experiment-2/scripts/compare_timed_track12.py"),
                    "--run-dir", str(run)], check=True)
    old = pd.read_csv(baseline / "results/track12_pointwise_comparison.csv")
    new = pd.read_csv(run / "results/track12_pointwise_comparison.csv")
    merged = old.merge(new, on=["species", "pair", "reference_step"], suffixes=("_old", "_new"), validate="one_to_one")
    if len(merged) != 13012 or len(old) != len(new):
        raise ValueError("Reference-grid samples differ")
    if not (merged.particle_id_old == merged.particle_id_new).all():
        raise ValueError("Particle identity mapping differs between the two comparison tables")
    delta_r = np.column_stack([merged[f"{c}_opalx_m_new"] - merged[f"{c}_opalx_m_old"] for c in ("x", "y", "s")])
    delta_p = np.column_stack([merged[f"p{c}_opalx_new"] - merged[f"p{c}_opalx_old"] for c in ("x", "y", "z")])
    merged["delta_r_m"] = np.linalg.norm(delta_r, axis=1)
    merged["delta_p_mc"] = np.linalg.norm(delta_p, axis=1)
    merged[["species", "pair", "reference_step", "particle_id_new", "delta_r_m", "delta_p_mc"]].to_csv(output / "trajectory_differences.csv", index=False)
    make_plot(merged, output)
    old_report = json.loads((baseline / "results/track12_comparison.json").read_text())
    new_report = json.loads((run / "results/track12_comparison.json").read_text())
    runtime = float(re.search(r"^real ([0-9.]+)$", (run / "opalx.log").read_text(), re.M)[1])
    old_runtime = float(re.search(r"^real ([0-9.]+)$", (baseline / "opalx.log").read_text(), re.M)[1])
    binary_hashes = (HERE / "build/binaries.sha256").read_text()
    executable_hash = re.search(r"^([0-9a-f]{64})\s+\S*/bin/opalx$", binary_hashes, re.M)[1]
    if executable_hash not in (run / "execution-fingerprint.txt").read_text():
        raise ValueError("Executed binary hash differs from the validated build")
    norms = paired_norms(baseline, run)
    summary = dict(validation=checks, baseline_validation=old_checks, input_hashes_match=True,
                   h5_ids_clocks_match=True,
                   raw_h5_bitwise_equal=all(a["canonical_h5_sha256"] == b["canonical_h5_sha256"]
                                           for a, b in zip(checks, old_checks)),
                   raw_h5_exact_numeric_equal=bool((stats.max_abs == 0).all()),
                   reference_samples=len(merged), **norms,
                   old_wall_seconds=old_runtime, new_wall_seconds=runtime,
                   old_cain=old_report["overall"], new_cain=new_report["overall"],
                   old_manufactured=old_report["opalx_vs_manufactured"],
                   new_manufactured=new_report["opalx_vs_manufactured"],
                   cutoff_patch_sha256=digest(HERE / "cutoff.patch"),
                   executable_sha256=executable_hash,
                   source_base_revision="e24dc92185198c5f65718e057400e6c6398dce56",
                   analysis_script_sha256=digest(Path(__file__)),
                   baseline=str(baseline), run=str(run),
                   acceptance="No tolerance inferred; report exact differences and unchanged-input evidence")
    if args.control:
        control = args.control.resolve()
        if not (control / "completed").is_file():
            raise ValueError("Old-binary repeat has not completed")
        for name, expected in preparation["input_sha256"].items():
            if digest(control / name) != expected:
                raise ValueError(f"Control input changed: {name}")
        control_checks = integrity(control, ip)
        repeatability = {"validation": control_checks,
                         "old_executable_sha256": "31e97e63d5208ef210c4188b717b48ccd739fba37874f38735e08754372d1307"}
        if repeatability["old_executable_sha256"] not in (control / "execution-fingerprint.txt").read_text():
            raise ValueError("Wrong old-binary control executable")
        for name, left, right in (("baseline_vs_old_repeat", baseline, control),
                                  ("old_repeat_vs_cutoff", control, run)):
            frame = comparer.compare(left, right)
            frame.to_csv(output / f"{name}.csv", index=False)
            repeatability[name] = paired_norms(left, right)
        repeatability["wall_seconds"] = float(re.search(r"^real ([0-9.]+)$", (control / "opalx.log").read_text(), re.M)[1])
        summary["old_binary_repeatability"] = repeatability
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(stats.to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
