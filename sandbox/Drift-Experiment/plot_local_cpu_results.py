#!/usr/bin/env python3
"""Plot BeamBeam witness results from H5/stat output."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUN_DIR = (
    ROOT
    / "sandbox"
    / "Drift-Experiment"
    / "one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs_selfsc_mesh24mm_nxy64"
)

ELECTRON_COLOR = "#b51f2e"
POSITRON_COLOR = "#1f5fb5"
LINE_STYLES = [
    "-",
    (0, (7, 5)),
    (0, (4, 4)),
    (0, (9, 6)),
    (0, (2, 3)),
    (0, (9, 4, 2, 4)),
]


def configure_matplotlib(output_dir: Path) -> None:
    cache_root = output_dir / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "figure.titlesize": 14,
            "legend.fontsize": 8,
            "lines.linewidth": 1.4,
            "savefig.dpi": 220,
        }
    )


def read_sdds_stat(path: Path) -> pd.DataFrame:
    text = path.read_text(encoding="utf-8")
    columns = re.findall(r"&column\s+name=([^,\n]+),", text)
    if not columns:
        raise ValueError(f"{path} has no SDDS column definitions")

    rows: list[list[str]] = []
    number = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
    for line in text.splitlines():
        values = line.split()
        if len(values) == len(columns) and number.match(values[0]):
            rows.append(values)

    data = pd.DataFrame(rows, columns=columns)
    for column in data.columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data


def sorted_steps(h5file: h5py.File) -> list[str]:
    return sorted(h5file.keys(), key=lambda key: int(key.split("#")[1]))


def reference_z(group: h5py.Group) -> float:
    if "RefPartR" in group.attrs:
        return float(group.attrs["RefPartR"][2])
    if "SPOS" in group.attrs:
        return float(group.attrs["SPOS"][0])
    return 0.0


def load_witness(path: Path, container: str, species: str, ip_s: float) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    with h5py.File(path, "r") as h5:
        for step_index, step_name in enumerate(sorted_steps(h5)):
            group = h5[step_name]
            time_ps = float(group.attrs["TIME"][0]) * 1.0e12
            ref_z = reference_z(group)
            ids = group["id"][:]
            ex = group["Ex"][:] if "Ex" in group else np.zeros_like(group["x"][:])
            ey = group["Ey"][:] if "Ey" in group else np.zeros_like(group["x"][:])
            ez = group["Ez"][:] if "Ez" in group else np.zeros_like(group["x"][:])
            for index, pid in enumerate(ids):
                pair = int(pid) + 1
                x = float(group["x"][index])
                y = float(group["y"][index])
                z_abs = ref_z + float(group["z"][index])
                rows.append(
                    {
                        "container": container,
                        "species": species,
                        "pair": pair,
                        "step": step_index,
                        "time_ps": time_ps,
                        "s_minus_ip_mm": 1.0e3 * (z_abs - ip_s),
                        "x_um": 1.0e6 * x,
                        "y_um": 1.0e6 * y,
                        "z_mm": 1.0e3 * z_abs,
                        "E_abs": float(np.sqrt(ex[index] ** 2 + ey[index] ** 2 + ez[index] ** 2)),
                        "Ex": float(ex[index]),
                        "Ey": float(ey[index]),
                        "Ez": float(ez[index]),
                    }
                )
    return pd.DataFrame(rows)


def plot_c0_envelope(c0: pd.DataFrame, output: Path, retire_ps: float, title_label: str) -> None:
    import matplotlib.pyplot as plt

    active = c0.loc[c0["numParticles"] > 0].copy()
    active["t_ps"] = active["t"] * 1.0e3

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.0), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.plot(active["t_ps"], active["rms_x"] * 1.0e3, label="rms_x", color="#1f77b4")
    ax.plot(active["t_ps"], active["rms_y"] * 1.0e3, label="rms_y", color="#ff7f0e")
    ax.plot(active["t_ps"], active["max_x"] * 1.0e3, "--", label="max_x", color="#1f77b4", alpha=0.8)
    ax.plot(active["t_ps"], active["max_y"] * 1.0e3, "--", label="max_y", color="#ff7f0e", alpha=0.8)
    ax.axhline(12.0, color="0.2", ls=":", lw=1.2, label="aperture half-width")
    ax.axvline(retire_ps, color="tab:green", ls="--", lw=1.1)
    ax.set_ylabel("transverse size [mm]")
    ax.set_title(f"{title_label}: c0 transverse expansion")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3)

    ax = axes[1]
    ax.plot(active["t_ps"], active["rms_s"] * 1.0e3, label="rms_s", color="#2ca02c")
    ax.plot(active["t_ps"], active["max_s"] * 1.0e3, "--", label="max_s", color="#2ca02c", alpha=0.8)
    ax.axvline(retire_ps, color="tab:green", ls="--", lw=1.1)
    ax.set_ylabel("longitudinal size [mm]")
    ax.set_title(f"{title_label}: c0 longitudinal expansion")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[2]
    ax.plot(active["t_ps"], active["energy"], color="#9467bd", label="c0 energy")
    ax2 = ax.twinx()
    ax2.step(c0["t"] * 1.0e3, c0["numParticles"], where="post", color="#555555", alpha=0.75, label="N")
    ax.axvline(retire_ps, color="tab:green", ls="--", lw=1.1)
    ax.set_xlabel("time [ps]")
    ax.set_ylabel("energy [MeV]")
    ax2.set_ylabel("macroparticles")
    ax.set_title(f"{title_label}: c0 energy and retirement")
    ax.grid(True, alpha=0.25)
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc="best")

    fig.savefig(output)
    plt.close(fig)


def plot_witness_trajectories(
    witness: pd.DataFrame, output: Path, retire_ps: float, title_label: str
) -> None:
    import matplotlib.pyplot as plt

    t0 = witness["time_ps"].min()
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True, constrained_layout=True)
    colors = {"c1": ELECTRON_COLOR, "c2": POSITRON_COLOR}
    labels = {"c1": "c1/e-", "c2": "c2/e+"}
    for container, group_c in witness.groupby("container", sort=False):
        for pair, group in group_c.groupby("pair", sort=True):
            style = LINE_STYLES[(int(pair) - 1) % len(LINE_STYLES)]
            label = f"{labels[container]} pair {pair}"
            axes[0].plot(
                group["time_ps"] - t0,
                group["x_um"],
                color=colors[container],
                ls=style,
                alpha=0.9,
                label=label,
            )
            axes[1].plot(
                group["time_ps"] - t0,
                group["y_um"],
                color=colors[container],
                ls=style,
                alpha=0.9,
            )
    for ax in axes:
        ax.axvline(retire_ps - t0, color="tab:green", ls="--", lw=1.1, label="c0 retire")
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("x [um]")
    axes[1].set_ylabel("y [um]")
    axes[1].set_xlabel("time since first witness dump [ps]")
    axes[0].set_title(f"{title_label}: witness x trajectories")
    axes[1].set_title(f"{title_label}: witness y trajectories")
    axes[0].legend(ncol=3, fontsize=7)
    fig.savefig(output)
    plt.close(fig)


def plot_witness_x_vs_s(witness: pd.DataFrame, output: Path, title_label: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    colors = {"c1": ELECTRON_COLOR, "c2": POSITRON_COLOR}
    labels = {"c1": "c1/e-", "c2": "c2/e+"}
    for container, group_c in witness.groupby("container", sort=False):
        for pair, group in group_c.groupby("pair", sort=True):
            style = LINE_STYLES[(int(pair) - 1) % len(LINE_STYLES)]
            ax.plot(
                group["s_minus_ip_mm"],
                group["x_um"],
                color=colors[container],
                ls=style,
                alpha=0.9,
                label=f"{labels[container]} pair {pair}",
            )
    ax.axhline(0.0, color="0.2", ls=(0, (2, 6)), lw=0.9)
    ax.set_xlabel("S - S_IP [mm]")
    ax.set_ylabel("x [um]")
    ax.set_title(f"{title_label}: Track12 witness x trajectories")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3, fontsize=7)
    fig.savefig(output)
    plt.close(fig)


def plot_witness_fields(witness: pd.DataFrame, output: Path, retire_ps: float, title_label: str) -> None:
    import matplotlib.pyplot as plt

    t0 = witness["time_ps"].min()
    fig, ax = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    colors = {"c1": ELECTRON_COLOR, "c2": POSITRON_COLOR}
    labels = {"c1": "c1/e-", "c2": "c2/e+"}
    for container, group_c in witness.groupby("container", sort=False):
        for pair, group in group_c.groupby("pair", sort=True):
            style = LINE_STYLES[(int(pair) - 1) % len(LINE_STYLES)]
            ax.semilogy(
                group["time_ps"] - t0,
                np.maximum(group["E_abs"], 1.0),
                color=colors[container],
                ls=style,
                alpha=0.65,
            )
    for container, group in witness.groupby("container", sort=False):
        summary = group.groupby("time_ps", as_index=False)["E_abs"].median()
        ax.semilogy(
            summary["time_ps"] - t0,
            np.maximum(summary["E_abs"], 1.0),
            color=colors[container],
            lw=2.4,
            label=f"{labels[container]} median",
        )
    ax.axvline(retire_ps - t0, color="tab:green", ls="--", lw=1.1, label="c0 retire")
    ax.set_xlabel("time since first witness dump [ps]")
    ax.set_ylabel("|E| sampled by witness [V/m]")
    ax.set_title(f"{title_label}: witness-sampled electric field")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(output)
    plt.close(fig)


def plot_c0_xy_snapshots(
    c0_h5: Path, output: Path, aperture_radius_m: float, sample: int, title_label: str
) -> None:
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(20260706)
    with h5py.File(c0_h5, "r") as h5:
        steps = sorted_steps(h5)
        selected = [steps[0], steps[-1]]
        snapshots = []
        for step in selected:
            group = h5[step]
            n = len(group["x"])
            indices = np.arange(n)
            if n > sample:
                indices = np.sort(rng.choice(indices, size=sample, replace=False))
            snapshots.append(
                {
                    "step": step,
                    "time_ps": float(group.attrs["TIME"][0]) * 1.0e12,
                    "x_mm": group["x"][indices] * 1.0e3,
                    "y_mm": group["y"][indices] * 1.0e3,
                    "rms_x_mm": float(group.attrs["RMSX"][0]) * 1.0e3,
                    "rms_y_mm": float(group.attrs["RMSX"][1]) * 1.0e3,
                }
            )

    theta = np.linspace(0.0, 2.0 * np.pi, 512)
    radius_mm = aperture_radius_m * 1.0e3
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.8), constrained_layout=True)
    for ax, snap in zip(axes, snapshots, strict=True):
        ax.plot(radius_mm * np.cos(theta), radius_mm * np.sin(theta), color="0.15", lw=1.1)
        ax.scatter(snap["x_mm"], snap["y_mm"], s=1.0, alpha=0.18, color="#1f77b4", rasterized=True)
        ellipse = plt.matplotlib.patches.Ellipse(
            (0.0, 0.0),
            2.0 * snap["rms_x_mm"],
            2.0 * snap["rms_y_mm"],
            fill=False,
            color="#d62728",
            lw=1.4,
            label="1 rms",
        )
        ax.add_patch(ellipse)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-radius_mm, radius_mm)
        ax.set_ylim(-radius_mm, radius_mm)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title(f"t = {snap['time_ps']:.1f} ps")
        ax.grid(True, alpha=0.18)
    axes[0].legend(loc="upper right")
    fig.suptitle(f"{title_label}: sampled c0 transverse distribution ({sample} particles per panel)")
    fig.savefig(output)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--ip-s", type=float, default=6.0e-3)
    parser.add_argument("--witness-t0-ps", type=float, default=4.0)
    parser.add_argument("--retire-ps", type=float, default=54.0)
    parser.add_argument("--aperture-radius-m", type=float, default=12.0e-3)
    parser.add_argument("--sample", type=int, default=30000)
    parser.add_argument("--title-label", default="CPU results")
    parser.add_argument("--output-prefix", default="cpu_results")
    args = parser.parse_args()

    run_dir = args.run_dir
    configure_matplotlib(run_dir)

    c0_stat = read_sdds_stat(run_dir / "spacecharge_drift_withness_c0.stat")
    c1 = load_witness(run_dir / "spacecharge_drift_withness_c1.h5", "c1", "electron", args.ip_s)
    c2 = load_witness(run_dir / "spacecharge_drift_withness_c2.h5", "c2", "positron", args.ip_s)
    witness = pd.concat([c1, c2], ignore_index=True)

    outputs = [
        run_dir / f"{args.output_prefix}_c0_envelope.png",
        run_dir / f"{args.output_prefix}_c0_xy_snapshots.png",
        run_dir / f"{args.output_prefix}_witness_trajectories_time.png",
        run_dir / f"{args.output_prefix}_track12_x_vs_s.png",
        run_dir / f"{args.output_prefix}_witness_fields.png",
    ]

    plot_c0_envelope(c0_stat, outputs[0], args.retire_ps, args.title_label)
    plot_c0_xy_snapshots(
        run_dir / "spacecharge_drift_withness_c0.h5",
        outputs[1],
        args.aperture_radius_m,
        args.sample,
        args.title_label,
    )
    plot_witness_trajectories(witness, outputs[2], args.retire_ps, args.title_label)
    plot_witness_x_vs_s(witness, outputs[3], args.title_label)
    plot_witness_fields(witness, outputs[4], args.retire_ps, args.title_label)

    for output in outputs:
        print(output.resolve())

    last_active = c0_stat.loc[c0_stat["numParticles"] > 0].iloc[-1]
    print(
        "last active c0: "
        f"t={last_active['t'] * 1.0e3:.3f} ps, "
        f"max_x={last_active['max_x'] * 1.0e3:.3f} mm, "
        f"max_y={last_active['max_y'] * 1.0e3:.3f} mm"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
