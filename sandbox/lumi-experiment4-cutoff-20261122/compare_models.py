#!/usr/bin/env python3
"""Compare existing source evolution and witness birth clocks; no new tracking.

Experiment4 primary snapshots come from the earlier two-rank CPU run because
the matched one-GPU run suppressed primary H5 output. Do not label these as
GPU measurements. Track12 source widths are computed by exact free drift of
the fixed input, justified by zero collective kicks and negligible divergence.
"""
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
C = 299792458.
IP = .17
T_IP = IP/C
TRACK = ROOT / "sandbox/lumi-cutoff-20261122/track12"
CPU = ROOT / "sandbox/experiment-4/results/cutoff-20261122-n2"
GPU = HERE / "cutoff"
OUT = HERE / "model-comparison"
STEM = "gamma_gamma_pairs-20mm-dynamic-xy-staged"


def scalar(g, attr):
    return float(np.asarray(g.attrs[attr]).item())


def snapshots(path):
    rows = []
    with h5py.File(path) as h:
        for g in h.values():
            if "x" not in g or len(g["x"]) == 0:
                continue
            if np.any(g.attrs["TaitBryantAngles"] != 0):
                raise ValueError("Nonzero rotation requires explicit transformation")
            r = np.column_stack([g[k][:] for k in ("x", "y", "z")]) + np.asarray(g.attrs["RefPartR"]).reshape(3)
            p = np.column_stack([g[k][:] for k in ("px", "py", "pz")])
            radius = np.linalg.norm(r[:, :2], axis=1)
            rows.append(dict(time_s=scalar(g, "TIME"), n=len(r),
                             sigma_x_m=np.std(r[:, 0]), sigma_y_m=np.std(r[:, 1]),
                             sigma_z_m=np.std(r[:, 2]), mean_s_m=r[:, 2].mean(),
                             median_r_m=np.median(radius), q90_r_m=np.quantile(radius, .9),
                             rms_px_mc=np.std(p[:, 0]), rms_py_mc=np.std(p[:, 1])))
    return pd.DataFrame(rows).sort_values("time_s")


def nearest(frame, time):
    return frame.iloc[(frame.time_s-time).abs().argmin()].to_dict()


def main():
    OUT.mkdir(exist_ok=True)
    primary = snapshots(CPU / f"{STEM}_c0.h5")
    primary.to_csv(OUT / "experiment4_primary_cpu.csv", index=False)
    at_ip = nearest(primary, T_IP)
    crossing = np.interp(IP, primary.mean_s_m, primary.time_s)
    fixed = np.loadtxt(TRACK / "input/primary_fixed.fromfile", skiprows=2)
    p = fixed[:, 3:6]
    v = C*p/np.sqrt(1+(p*p).sum(axis=1))[:, None]
    flight = np.linspace(0, 1501*1.8e-6/C, 30)
    rigid = pd.DataFrame([dict(time_s=t, sigma_x_m=np.std(fixed[:, 0]+v[:, 0]*t),
                               sigma_y_m=np.std(fixed[:, 1]+v[:, 1]*t)) for t in flight])
    rigid.to_csv(OUT / "track12_predicted_primary_free_drift.csv", index=False)
    pairs = np.loadtxt(TRACK / "input/track12_positrons.emittedfromfile", skiprows=2)
    births = pd.DataFrame(dict(pair=np.arange(1, 7), birth_relative_nominal_overlap_ps=pairs[:, 6]*1e12,
                               z_at_birth_relative_ip_mm=pairs[:, 2]*1e3))
    births.to_csv(OUT / "track12_births.csv", index=False)
    summaries = dict(primary_cpu_diagnostic=str(CPU / f"{STEM}_c0.h5"),
                     track12_primary_input=str(TRACK / "input/primary_fixed.fromfile"),
                     nominal_overlap_time_s=T_IP,
                     experiment4_measured_cpu_centroid_crossing_s=float(crossing),
                     experiment4_primary_near_nominal_overlap=at_ip,
                     experiment4_primary_near_witness_injection=nearest(primary, T_IP-16.747e-12),
                     track12_primary_input_sigma_x_m=float(rigid.sigma_x_m.iloc[0]),
                     track12_primary_input_sigma_y_m=float(rigid.sigma_y_m.iloc[0]),
                     track12_sigma_x_change_over_run_m=float(rigid.sigma_x_m.iloc[-1]-rigid.sigma_x_m.iloc[0]),
                     track12_initial_rms_px_mc=float(np.std(p[:, 0])),
                     track12_initial_rms_py_mc=float(np.std(p[:, 1])),
                     experiment4_input_nominal_sigma_xy_m=float(np.sqrt(5.2e-9*.727e-3)),
                     experiment4_input_nominal_sigma_xprime_rad=float(np.sqrt(5.2e-9/.727e-3)),
                     experiment4_unforced_emittance_drift_sigma_at_ip_m=float(np.sqrt(5.2e-9*.727e-3+IP**2*5.2e-9/.727e-3)),
                     experiment4_injection_before_nominal_overlap_ps=-16.747,
                     experiment4_nominal_distance_per_centroid_at_injection_m=16.747e-12*C,
                     primary_sigma_z_over_c_ps=.6e-3/C*1e12,
                     track12_births=births.to_dict("records"))
    for c, species in ((1, "electron"), (2, "positron")):
        frame = snapshots(GPU / f"{STEM}_c{c}.h5")
        frame.to_csv(OUT / f"experiment4_gpu_{species}_envelope.csv", index=False)
        initial = np.loadtxt(GPU / f"input/gamma_gamma_{species}s.fromfile", skiprows=2)
        pp = initial[:, 3:6]
        vv = C*pp/np.sqrt(1+(pp*pp).sum(axis=1))[:, None]
        summaries[f"experiment4_{species}"] = dict(
            first_saved=frame.iloc[0].to_dict(), at_nominal_overlap=nearest(frame, T_IP),
            input_median_vtransverse_over_c=float(np.median(np.linalg.norm(vv[:, :2], axis=1))/C),
            unforced_median_transverse_travel_before_overlap_m=float(np.median(np.linalg.norm(vv[:, :2], axis=1))*16.747e-12))
    summaries["primary_width_ratio_x"] = at_ip["sigma_x_m"]/rigid.sigma_x_m.iloc[0]
    summaries["scope"] = "CPU2 Experiment4 source snapshots; GPU1 witnesses; Track12 source predicted by free drift. Not a matched intervention or measured GPU source comparison."
    (OUT / "summary.json").write_text(json.dumps(summaries, indent=2)+"\n")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    pre_and_overlap = primary.time_s <= T_IP+40e-12
    for component, style in (("x", "-"), ("y", "--")):
        axes[0].plot((primary.time_s[pre_and_overlap]-T_IP)*1e12, primary[f"sigma_{component}_m"][pre_and_overlap]*1e6,
                     style, color="#b51f2e", label=f"Experiment 4 CPU: σ{component}")
        close = abs(primary.time_s-T_IP) <= 20e-12
        axes[1].plot((primary.time_s[close]-T_IP)*1e12, primary[f"sigma_{component}_m"][close]*1e6,
                     style, color="#b51f2e")
    axes[0].set(title="Evolving primary: stored CPU snapshots", xlabel="Time from nominal overlap [ps]", ylabel="Primary RMS size [µm]", xlim=(-T_IP*1e12, 40), ylim=(0, 500))
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].hlines(rigid.sigma_x_m.iloc[0]*1e6, -20, 20, color="#1f5fb5", label="Track-12: free-drift prediction")
    axes[1].set(title="Source size near overlap", xlabel="Time from nominal overlap [ps]", ylabel="Primary RMS size [µm]", yscale="log")
    axes[1].legend(frameon=False, fontsize=8, loc="center right")
    axes[2].axvspan(-3*.6e-3/C*1e12, 3*.6e-3/C*1e12, alpha=.12, color="gray", label="±3 nominal σz/c")
    axes[2].scatter(births.birth_relative_nominal_overlap_ps, np.arange(1, 7), color="#1f5fb5", label="Track-12 births")
    axes[2].axvline(-16.747, color="#b51f2e", label="Experiment 4: all births")
    axes[2].axvline(0, color="k", ls=":", lw=1)
    axes[2].set(title="Witness injection timeline", xlabel="Time from nominal overlap [ps]", ylabel="Track-12 pair", xlim=(-20, 8))
    axes[2].legend(frameon=False, fontsize=8, loc="upper left")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(OUT / "primary_evolution_and_births.png", dpi=220, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
