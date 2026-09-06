#!/usr/bin/env python3
"""Run one-rank DT convergence against the OPAL 2022.1 RK4 reference.

Usage: ~/.venv-h6/bin/python compare_coasting72.py --old-orbit /path/to/cyclotron1-trackOrbit.dat
All runs use fresh directories under --output; the original inputs are preserved.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_run(path):
    rows = []
    with h5py.File(path) as file:
        for group in file.values():
            if "RefPartR" not in group.attrs:
                continue
            row = {"t": float(group.attrs["TIME"][0]), "s": float(group.attrs["SPOS"][0]),
                   "energy": float(group.attrs["ENERGY"][0])}
            for key, value in zip(("x", "y", "z"), group.attrs["RefPartR"]):
                row[key] = value
            for key, value in zip(("px", "py", "pz"), group.attrs["RefPartP"]):
                row[key] = value
            row["particle_offset"] = float(np.linalg.norm([group[k][0] for k in ("x", "y", "z")]))
            rows.append(row)
    return pd.DataFrame(rows).sort_values("t").reset_index(drop=True)


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-orbit", type=Path, default=root/"reference72/old-rk4-orbit.dat")
    parser.add_argument("--lf2-orbit", type=Path, default=root/"reference72/old-lf2-orbit.dat")
    parser.add_argument("--executable", type=Path, default=root.parents[1]/"omp-build/src/opalx")
    parser.add_argument("--output", type=Path, default=root/"results72")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    old = np.loadtxt(args.old_orbit, usecols=range(1, 7))
    old_r = old[:, [0, 4, 2]] * [1, -1, 1]
    old_p = old[:, [1, 5, 3]] * [1, -1, 1]
    template = (root/"opalx/coasting72.in").read_text()
    # This is the orbit-only regression; the separate map case is diagnostic.
    template = template.replace('ENABLELINEARTRANSFERMAPS=TRUE', 'ENABLELINEARTRANSFERMAPS=FALSE')
    dt = 1.64527805199081e-10
    metrics = []
    runs = {}
    for factor, ranks in [(1, 1), (2, 1), (4, 1)]:
        name = f"dt{factor}_rank{ranks}"
        work = args.output/name
        work.mkdir(exist_ok=False)
        text = template.replace('../opal/bfield.dat', str(root/"opal/bfield.dat"))
        text = text.replace('proton72.dat', str(root/"opalx/proton72.dat"))
        text = re.sub(r"DT=1.64527805199081e-10", f"DT={dt/factor:.17g}", text)
        text = text.replace("MAXSTEPS=2000", f"MAXSTEPS={2000*factor}")
        (work/"coasting72.in").write_text(text)
        command = [str(args.executable.resolve()), "--info", "2", "coasting72.in"]
        env = dict(os.environ, OMP_NUM_THREADS="2", OMP_PROC_BIND="false")
        with (work/"run.log").open("w") as log:
            subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        log = (work/"run.log").read_text()
        if "completed directed turn 1" not in log:
            raise RuntimeError(f"{name}: no directed turn completion")
        data = read_run(work/"coasting72.h5")
        runs[name] = data
        # H5 Step#0 is the first endpoint, while old orbit row 0 is the launch.
        sampled = data.iloc[factor-1:(len(old)-1)*factor:factor]
        r_error = np.linalg.norm(sampled[["x", "y", "z"]].to_numpy()-old_r[1:], axis=1)
        p_error = np.linalg.norm(sampled[["px", "py", "pz"]].to_numpy()-old_p[1:], axis=1)
        metrics.append(dict(run=name, dt_s=dt/factor, ranks=ranks,
            max_position_error_m=float(r_error.max()), max_momentum_error=float(p_error.max()),
            energy_drift_MeV=float(np.ptp(data.energy)), max_particle_offset_m=float(data.particle_offset.max()),
            end_time_s=float(data.t.iloc[-1]), end_path_m=float(data.s.iloc[-1]), steps=len(data)))
        data.to_csv(work/"reference.csv", index=False)
        pd.DataFrame({"t_s":sampled.t.to_numpy(), "position_error_m":r_error}).to_csv(work/"errors.csv", index=False)
        print(json.dumps(metrics[-1]), flush=True)
    result = pd.DataFrame(metrics)
    result.to_csv(args.output/"comparison.csv", index=False)
    lf2 = np.loadtxt(args.lf2_orbit, usecols=range(1, 7))
    baseline = runs["dt1_rank1"]
    n = min(len(baseline), len(lf2)-1)
    lf2_r = lf2[1:n+1, [0,4,2]] * [1,-1,1]
    lf2_p = lf2[1:n+1, [1,5,3]] * [1,-1,1]
    matched = dict(
        integrator="old OPAL LF2 versus OPALX Boris", ranks=1, samples=n,
        max_position_error_m=float(np.linalg.norm(baseline[["x","y","z"]].to_numpy()[:n]-lf2_r,axis=1).max()),
        max_momentum_error=float(np.linalg.norm(baseline[["px","py","pz"]].to_numpy()[:n]-lf2_p,axis=1).max()))
    (args.output/"matched-integrator.json").write_text(json.dumps(matched,indent=2)+"\n")
    print(json.dumps(matched),flush=True)
    # Legacy orbit ASCII has eight digits after the decimal, limiting position
    # comparison to several nanometres at metre radii. Allow 20 nm for that
    # rounding and the 9-digit initial momentum retained in the launch fixture.
    assert matched["max_position_error_m"] < 2e-8
    assert matched["max_momentum_error"] < 3e-9
    assert result.energy_drift_MeV.max() < 1e-9
    assert result.max_particle_offset_m.max() < 1e-9
    assert metrics[2]["max_position_error_m"] < metrics[0]["max_position_error_m"]/8
    plt.rcParams.update({"font.size":11, "axes.spines.top":False, "axes.spines.right":False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), layout="constrained")
    axes[0].plot(old_r[:,0], old_r[:,2], color="black", lw=1.5, label="OPAL 2022.1 RK4")
    fine = runs["dt4_rank1"]
    axes[0].plot(fine.x, fine.z, color="#0072B2", lw=1, ls="--", label="OPALX, DT/4")
    axes[0].set(xlabel="X [m]", ylabel="Z [m]", title="72 MeV coasting proton", aspect="equal")
    axes[0].legend(frameon=False)
    for factor in (1,2,4):
        errors = pd.read_csv(args.output/f"dt{factor}_rank1/errors.csv")
        axes[1].plot(errors.t_s*1e9, errors.position_error_m*1e6, label=f"DT/{factor}")
    axes[1].set(xlabel="Time [ns]", ylabel="Position difference [µm]", title="Difference from old OPAL")
    axes[1].legend(frameon=False)
    fig.savefig(args.output/"coasting72.png", dpi=220)
    fig.savefig(args.output/"coasting72.pdf")


if __name__ == "__main__":
    main()
