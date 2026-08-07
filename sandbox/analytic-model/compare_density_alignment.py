#!/usr/bin/env python3
"""Compare OPALX C0 density plots against MSOL plot-time field maps."""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "figures" / "opalx-msol-density-alignment"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def configure_matplotlib() -> None:
    cache_dir = SCRIPT_DIR / "atic" / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def normalize_shape(values: np.ndarray) -> np.ndarray:
    out = values.astype(float).copy()
    positive = out[np.isfinite(out) & (out > 0.0)]
    if positive.size == 0:
        return out
    out /= float(np.nanmax(positive))
    return out


def msol_on_xy(opal_frame: dict, msol_config: dict, time_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    make_plots = load_module("make_plots_for_density", SCRIPT_DIR / "make_plots.py")
    track12, sources = make_plots.make_msol_sources(msol_config)
    x_edges = opal_frame["xy_edges"]
    y_edges = opal_frame["xy_edges"]
    x_um = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_um = 0.5 * (y_edges[:-1] + y_edges[1:])
    values = np.empty((len(y_um), len(x_um)), dtype=float)
    for iy, y_value in enumerate(y_um):
        for ix, x_value in enumerate(x_um):
            position = np.array([x_value * 1.0e-6, y_value * 1.0e-6, 0.0])
            e_field, _ = track12.anisotropic_total_lab_fields(position, time_s, sources)
            values[iy, ix] = float(np.dot(e_field, e_field))
    return x_edges, y_edges, values


def msol_on_zr(opal_frame: dict, msol_config: dict, time_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    make_plots = load_module("make_plots_for_density", SCRIPT_DIR / "make_plots.py")
    track12, sources = make_plots.make_msol_sources(msol_config)
    z_edges = opal_frame["z_edges"]
    r_edges = opal_frame["r_edges"]
    z_um = 0.5 * (z_edges[:-1] + z_edges[1:])
    r_um = 0.5 * (r_edges[:-1] + r_edges[1:])
    values = np.empty((len(r_um), len(z_um)), dtype=float)
    for ir, r_value in enumerate(r_um):
        for iz, z_value in enumerate(z_um):
            position = np.array([r_value * 1.0e-6, 0.0, z_value * 1.0e-6])
            e_field, _ = track12.anisotropic_total_lab_fields(position, time_s, sources)
            values[ir, iz] = float(np.dot(e_field, e_field))
    return z_edges, r_edges, values


def msol_witness_overlay(trajectory: pd.DataFrame, step: int, projection: str) -> dict[str, np.ndarray]:
    frame = trajectory.loc[trajectory["step"] == step].copy()
    output: dict[str, np.ndarray] = {}
    for species, group in frame.groupby("species", sort=True):
        group = group.sort_values("pair")
        if projection == "zr":
            x = group["s_mm"].to_numpy(dtype=float) * 1000.0
            y = 1.0e6 * np.sqrt(group["x"].to_numpy(dtype=float) ** 2 + group["y"].to_numpy(dtype=float) ** 2)
        elif projection == "xy":
            x = 1.0e6 * group["x"].to_numpy(dtype=float)
            y = 1.0e6 * group["y"].to_numpy(dtype=float)
        else:
            raise ValueError(f"unknown projection {projection!r}")
        pair = group["pair"].to_numpy(dtype=float)
        output[species] = np.column_stack((pair, x, y))
    return output


def plot_three_panel(
    output: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    edges_x: np.ndarray,
    edges_y: np.ndarray,
    panels: list[tuple[str, np.ndarray, float]],
    overlays: list[dict[str, np.ndarray] | None] | None = None,
) -> None:
    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "savefig.dpi": 220,
        }
    )
    normalized = [normalize_shape(values) for _, values, _ in panels]
    positive = np.concatenate([values[np.isfinite(values) & (values > 0.0)] for values in normalized])
    vmin = max(float(np.nanmin(positive)), 1.0e-8)
    vmax = 1.0

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.6), sharex=True, sharey=True)
    mesh = None
    if overlays is None:
        overlays = [None] * len(panels)
    styles = {
        "electron": {"marker": "o", "edgecolors": "#c00000", "facecolors": "none", "label": "e-"},
        "positron": {"marker": "s", "edgecolors": "#0017b8", "facecolors": "none", "label": "e+"},
    }
    for ax, (label, _raw, time_ps), values, overlay in zip(axes, panels, normalized, overlays, strict=True):
        mesh = ax.pcolormesh(
            edges_x,
            edges_y,
            np.ma.masked_invalid(values),
            cmap="magma",
            norm=LogNorm(vmin=vmin, vmax=vmax),
            shading="auto",
        )
        if overlay:
            for species, points in overlay.items():
                style = styles[species]
                ax.scatter(
                    points[:, 1],
                    points[:, 2],
                    s=42 if species == "electron" else 30,
                    marker=style["marker"],
                    facecolors=style["facecolors"],
                    edgecolors=style["edgecolors"],
                    linewidths=1.1,
                    label=style["label"],
                    zorder=5,
                )
                for pair, x_value, y_value in points:
                    ax.annotate(
                        str(int(pair)),
                        (x_value, y_value),
                        xytext=(3, 3 if species == "electron" else -8),
                        textcoords="offset points",
                        fontsize=7,
                        color=style["edgecolors"],
                        zorder=6,
                    )
        ax.set_title(f"{label}\nt = {time_ps:.4g} ps")
        ax.set_xlabel(xlabel)
        ax.set_xlim(float(edges_x[0]), float(edges_x[-1]))
        ax.set_ylim(float(edges_y[0]), float(edges_y[-1]))
        ax.tick_params(direction="out", top=True, right=True)
    axes[0].set_ylabel(ylabel)
    fig.suptitle(title, y=0.98)
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes, fraction=0.025, pad=0.015)
        cbar.set_label(r"normalized $\mathbf{E}\cdot\mathbf{E}$")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles[:2], labels[:2], loc="upper right", bbox_to_anchor=(0.915, 0.94), fontsize=8)
    fig.text(
        0.5,
        0.01,
        f"directory: {Path.cwd().resolve()}    script: {Path(__file__).resolve()}",
        ha="center",
        va="bottom",
        fontsize=6,
        color="0.35",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.06, right=0.92, bottom=0.16, top=0.82, wspace=0.08)
    fig.savefig(output)
    plt.close(fig)


def main() -> int:
    opalc0 = load_module("opalx_c0_contours_for_density", SCRIPT_DIR / "opalx_c0_contours.py")
    make_plots = load_module("make_plots_for_density", SCRIPT_DIR / "make_plots.py")
    opal_config = opalc0.load_config(SCRIPT_DIR / "opalx_c0_contours.json")
    msol_config = make_plots.load_config(SCRIPT_DIR / "msol_matched_charge.json")
    trajectory = pd.read_csv(msol_config["trajectory_csv"])
    step_times = trajectory.groupby("step")["time_s"].first()

    with h5py.File(opal_config["h5_file"], "r") as handle:
        opal_frame = opalc0.load_step_grid(handle["Step#0"], opal_config)
    opal_overlay = opalc0.load_witness_frame(opal_config, 0)

    opal_time_ps = opal_frame["time_s"] * 1.0e12
    msol0_time_s = float(step_times.loc[0])
    msol500_time_s = float(step_times.loc[500])

    xy_edges, _, msol_xy0 = msol_on_xy(opal_frame, msol_config, msol0_time_s)
    _, _, msol_xy500 = msol_on_xy(opal_frame, msol_config, msol500_time_s)
    plot_three_panel(
        OUTPUT_DIR / "density-alignment-xy-opalx0-msol0-msol500.png",
        r"Transverse density-shape comparison: $\mathbf{E}\cdot\mathbf{E}$",
        r"x [$\mu$m]",
        r"y [$\mu$m]",
        xy_edges,
        xy_edges,
        [
            ("OPALX C0 step 0", opal_frame["xy_grid"], opal_time_ps),
            ("MSOL step 0", msol_xy0, msol0_time_s * 1.0e12),
            ("MSOL step 500", msol_xy500, msol500_time_s * 1.0e12),
        ],
        [
            None,
            msol_witness_overlay(trajectory, 0, "xy"),
            msol_witness_overlay(trajectory, 500, "xy"),
        ],
    )

    z_edges, r_edges, msol_zr0 = msol_on_zr(opal_frame, msol_config, msol0_time_s)
    _, _, msol_zr500 = msol_on_zr(opal_frame, msol_config, msol500_time_s)
    plot_three_panel(
        OUTPUT_DIR / "density-alignment-zr-opalx0-msol0-msol500.png",
        r"Longitudinal-radial density-shape comparison: $\mathbf{E}\cdot\mathbf{E}$",
        r"z - IP [$\mu$m]",
        r"r [$\mu$m]",
        z_edges,
        r_edges,
        [
            ("OPALX C0 step 0", opal_frame["zr_grid"], opal_time_ps),
            ("MSOL step 0", msol_zr0, msol0_time_s * 1.0e12),
            ("MSOL step 500", msol_zr500, msol500_time_s * 1.0e12),
        ],
        [
            opal_overlay,
            msol_witness_overlay(trajectory, 0, "zr"),
            msol_witness_overlay(trajectory, 500, "zr"),
        ],
    )

    print(f"OPALX step 0 time: {opal_time_ps:.6g} ps")
    print(f"MSOL step 0 time: {msol0_time_s * 1.0e12:.6g} ps")
    print(f"MSOL step 500 time: {msol500_time_s * 1.0e12:.6g} ps")
    print(f"Wrote {OUTPUT_DIR / 'density-alignment-xy-opalx0-msol0-msol500.png'}")
    print(f"Wrote {OUTPUT_DIR / 'density-alignment-zr-opalx0-msol0-msol500.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
