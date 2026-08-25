#!/usr/bin/env python3
"""Generate PNG figures for the CAIN-only manufactured analytic solution."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CAIN_CSV = ROOT / "sandbox" / "track12particles" / "outputs" / "track12_trajectory.csv"
DEFAULT_DATA_DIR = SCRIPT_DIR / "data"
DEFAULT_TRAJECTORY = DEFAULT_DATA_DIR / "cain-fit-half-charge-analytic-trajectory.csv"
DEFAULT_FIGURE_DIR = SCRIPT_DIR / "figures"
DEFAULT_CONFIG = SCRIPT_DIR / "manufactured_sol.json"
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"

ELECTRON_COLOR = "#c00000"
POSITRON_COLOR = "#0017b8"
LINE_STYLES = {
    1: "-",
    2: (0, (8, 5)),
    3: (0, (5, 5)),
    4: (0, (9, 6)),
    5: (0, (2, 3)),
    6: (0, (10, 4, 2, 4)),
}


def resolve_config_path(config_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return config_dir / path


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config_dir = path.resolve().parent
    return {
        "config_path": path.resolve(),
        "cain_csv": resolve_config_path(config_dir, config.get("cain_csv")) or DEFAULT_CAIN_CSV,
        "data_dir": resolve_config_path(config_dir, config.get("data_dir")) or DEFAULT_DATA_DIR,
        "figure_dir": resolve_config_path(config_dir, config.get("figure_dir")) or DEFAULT_FIGURE_DIR,
        "trajectory_csv": resolve_config_path(config_dir, config.get("trajectory_csv")) or DEFAULT_TRAJECTORY,
        "field_summary_csv": resolve_config_path(config_dir, config.get("field_summary_csv")),
        "trajectory_plot": resolve_config_path(config_dir, config.get("trajectory_plot")),
        "field_plot": resolve_config_path(config_dir, config.get("field_plot")),
        "field_components_step0_plot": resolve_config_path(config_dir, config.get("field_components_step0_plot")),
        "field_xy_step0_contour_plot": resolve_config_path(config_dir, config.get("field_xy_step0_contour_plot")),
        "field_zr_step0_contour_plot": resolve_config_path(config_dir, config.get("field_zr_step0_contour_plot")),
        "field_xy_sequence_dir": resolve_config_path(config_dir, config.get("field_xy_sequence_dir")),
        "field_zr_sequence_dir": resolve_config_path(config_dir, config.get("field_zr_sequence_dir")),
        "field_contour_steps": int(config.get("field_contour_steps", 0)),
        "field_contour_stride": int(config.get("field_contour_stride", 1)),
        "charge_scale": float(config.get("charge_scale", 0.5)),
        "source_selection": config.get("source_selection", "oncoming"),
        "source_time_offset_ps": float(config.get("source_time_offset_ps", 0.0)),
        "source_sigma_z_m": float(config.get("source_sigma_z_m", 0.6e-3)),
        "mc_source_particles": int(config.get("mc_source_particles", 400000)),
        "field_xy_extent_um": float(config.get("field_xy_extent_um", 8.0)),
        "field_xy_points": int(config.get("field_xy_points", 121)),
        "field_z_extent_um": float(config.get("field_z_extent_um", 1000.0)),
        "field_r_extent_um": float(config.get("field_r_extent_um", 8.0)),
        "field_zr_points": int(config.get("field_zr_points", 121)),
    }


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_msol_sources(config: dict):
    track12 = load_track12_module()
    source_charge_C = track12.DEFAULT_SOURCE_CHARGE_C * config["charge_scale"]
    sources = track12.make_anisotropic_sources(
        track12.DEFAULT_SOURCE_KINETIC_MEV,
        source_charge_C,
        track12.SIGMA_XY_M,
        track12.SIGMA_XY_M,
        config["source_sigma_z_m"],
        config["source_selection"],
        config["source_time_offset_ps"] * 1.0e-12,
        config["mc_source_particles"],
        track12.DEFAULT_MC_SEED,
    )
    return track12, sources


def evaluate_xy_edot(config: dict, time_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    track12, sources = make_msol_sources(config)
    extent_um = config["field_xy_extent_um"]
    points = config["field_xy_points"]
    if points < 11:
        raise ValueError("field_xy_points must be at least 11")

    x_um = np.linspace(-extent_um, extent_um, points)
    y_um = np.linspace(-extent_um, extent_um, points)
    e_dot_e = np.empty((points, points), dtype=float)
    for iy, y_value_um in enumerate(y_um):
        for ix, x_value_um in enumerate(x_um):
            position = np.array([x_value_um * 1.0e-6, y_value_um * 1.0e-6, 0.0], dtype=float)
            e_field, _b_field = track12.anisotropic_total_lab_fields(position, time_s, sources)
            e_dot_e[iy, ix] = float(np.dot(e_field, e_field))
    return x_um, y_um, e_dot_e


def evaluate_zr_edot(config: dict, time_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    track12, sources = make_msol_sources(config)
    z_extent_um = config["field_z_extent_um"]
    r_extent_um = config["field_r_extent_um"]
    points = config["field_zr_points"]
    if points < 11:
        raise ValueError("field_zr_points must be at least 11")

    z_um = np.linspace(-z_extent_um, z_extent_um, points)
    r_um = np.linspace(0.0, r_extent_um, points)
    e_dot_e = np.empty((points, points), dtype=float)
    for ir, r_value_um in enumerate(r_um):
        for iz, z_value_um in enumerate(z_um):
            position = np.array([r_value_um * 1.0e-6, 0.0, z_value_um * 1.0e-6], dtype=float)
            e_field, _b_field = track12.anisotropic_total_lab_fields(position, time_s, sources)
            e_dot_e[ir, iz] = float(np.dot(e_field, e_field))
    return z_um, r_um, e_dot_e


def edot_levels(frames: list[np.ndarray]) -> np.ndarray:
    positive_parts = [frame[frame > 0.0] for frame in frames]
    positive = np.concatenate([part for part in positive_parts if part.size])
    if positive.size == 0:
        return np.asarray([])
    return np.geomspace(float(positive.min()), float(positive.max()), 42)


def configure_matplotlib(output_dir: Path) -> None:
    cache_dir = SCRIPT_DIR / "atic" / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)


def add_footer(fig, directory: Path, config_path: Path) -> None:
    wrap_width = max(80, int(fig.get_figwidth() * 14))
    footer_text = "\n".join(
        [
            textwrap.fill(f"directory: {directory}", width=wrap_width),
            textwrap.fill(f"config: {config_path}", width=wrap_width),
        ]
    )
    fig.text(0.5, 0.01, footer_text, ha="center", va="bottom", fontsize=5.8, color="0.35")


def plot_trajectory_comparison(
    reference: pd.DataFrame,
    manufactured: pd.DataFrame,
    output: Path,
    directory: Path,
    config_path: Path,
    source_selection: str,
    charge_scale: float,
) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 18,
            "axes.titlesize": 16,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "savefig.dpi": 220,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.1), sharex=True, sharey=True)
    panels = (
        (axes[0], reference, "CAIN reference"),
        (axes[1], manufactured, "Manufactured model"),
    )
    x_min = min(reference["s_mm"].min(), manufactured["s_mm"].min()) - 0.05
    x_max = max(reference["s_mm"].max(), manufactured["s_mm"].max()) + 0.05
    for ax, frame, title in panels:
        for (_name, pair, kind), group in frame.groupby(["name", "pair", "kind"], sort=False):
            color = ELECTRON_COLOR if kind == 2 else POSITRON_COLOR
            group = group.sort_values("step")
            ax.plot(
                group["s_mm"],
                group["x_um"],
                color=color,
                linewidth=0.9,
                linestyle=LINE_STYLES.get(int(pair), "-"),
            )
            initial = group.iloc[0]
            ax.plot(
                initial["s_mm"],
                initial["x_um"],
                marker="o",
                markersize=2.8,
                markerfacecolor="white",
                markeredgecolor=color,
                markeredgewidth=0.7,
                linestyle="None",
            )
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle=(0, (2, 6)))
        ax.set_title(title)
        ax.set_xlabel("S (mm)")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-9.0, 20.0)
        ax.tick_params(direction="out", top=True, right=True, which="both")
        ax.minorticks_on()

    axes[0].set_ylabel(r"x ($\mu$m)")
    fig.suptitle(f"12 Test Particle Trajectories ({source_selection}, charge scale {charge_scale:g})", y=0.98)
    fig.tight_layout(rect=(0.0, 0.12, 1.0, 0.94))
    add_footer(fig, directory, config_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_field_vs_time(
    trajectory: pd.DataFrame,
    output: Path,
    summary_path: Path,
    directory: Path,
    config_path: Path,
) -> Path:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    required = {"species", "pair", "time_s", "E_abs_V_per_m"}
    missing = sorted(required.difference(trajectory.columns))
    if missing:
        raise ValueError(f"trajectory is missing columns: {', '.join(missing)}")

    data = trajectory.copy()
    data["time_ps"] = data["time_s"] * 1.0e12

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 9,
            "lines.linewidth": 1.15,
            "savefig.dpi": 220,
        }
    )

    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    for (species, pair), group in data.groupby(["species", "pair"], sort=True):
        color = ELECTRON_COLOR if species == "electron" else POSITRON_COLOR
        group = group.sort_values("time_s")
        ax.semilogy(
            group["time_ps"],
            group["E_abs_V_per_m"],
            color=color,
            ls=LINE_STYLES.get(int(pair), "-"),
            alpha=0.92,
        )

    ax.axvline(0.0, color="black", lw=0.9, ls=(0, (2, 6)))
    ax.set_xlabel("time [ps]")
    ax.set_ylabel(r"$|E|$ at particle position [V/m]")
    ax.set_title("Analytic half-charge oncoming source: sampled field")
    ax.grid(True, which="both", alpha=0.25)
    ax.minorticks_on()
    species_handles = [
        plt.Line2D([0], [0], color=ELECTRON_COLOR, lw=2.0, label="e-"),
        plt.Line2D([0], [0], color=POSITRON_COLOR, lw=2.0, label="e+"),
    ]
    pair_handles = [
        plt.Line2D([0], [0], color="0.25", lw=1.2, ls=LINE_STYLES[pair], label=f"pair {pair}")
        for pair in range(1, 7)
    ]
    first = ax.legend(handles=species_handles, loc="upper left", title="species")
    ax.add_artist(first)
    ax.legend(handles=pair_handles, loc="upper right", title="timing")
    fig.tight_layout(rect=(0.0, 0.15, 1.0, 0.95))
    add_footer(fig, directory, config_path)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)

    summary = (
        data.groupby(["species", "pair"], as_index=False)
        .agg(
            samples=("E_abs_V_per_m", "size"),
            time_min_ps=("time_ps", "min"),
            time_max_ps=("time_ps", "max"),
            E_min_V_per_m=("E_abs_V_per_m", "min"),
            E_max_V_per_m=("E_abs_V_per_m", "max"),
        )
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    return summary_path


def plot_step0_field_components(
    trajectory: pd.DataFrame,
    output: Path,
    directory: Path,
    config_path: Path,
) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    required = {"species", "pair", "step", "s_mm", "Ex_V_per_m", "Ey_V_per_m", "Ez_V_per_m"}
    missing = sorted(required.difference(trajectory.columns))
    if missing:
        raise ValueError(f"trajectory is missing columns: {', '.join(missing)}")

    data = trajectory.loc[trajectory["step"] == 0].copy()
    if data.empty:
        raise ValueError("trajectory has no step-0 rows")
    data = data.sort_values(["species", "pair"])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "savefig.dpi": 220,
        }
    )

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.4), sharex=True)
    components = (
        ("Ex_V_per_m", r"$E_x$ [V/m]"),
        ("Ey_V_per_m", r"$E_y$ [V/m]"),
        ("Ez_V_per_m", r"$E_z$ [V/m]"),
    )
    markers = {"electron": "o", "positron": "s"}
    colors = {"electron": ELECTRON_COLOR, "positron": POSITRON_COLOR}
    for ax, (column, ylabel) in zip(axes, components, strict=True):
        for species, group in data.groupby("species", sort=True):
            ax.plot(
                group["s_mm"],
                group[column],
                marker=markers.get(species, "o"),
                markersize=7.0 if species == "electron" else 4.5,
                linestyle="None",
                color=colors.get(species, "0.2"),
                markerfacecolor="none" if species == "electron" else colors.get(species, "0.2"),
                markeredgewidth=1.1 if species == "electron" else 0.8,
                label="e-" if species == "electron" else "e+",
                zorder=3 if species == "electron" else 4,
            )
            for _, row in group.iterrows():
                ax.annotate(
                    str(int(row["pair"])),
                    (row["s_mm"], row[column]),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=7,
                    color=colors.get(species, "0.2"),
                )
        ax.axhline(0.0, color="black", lw=0.8, ls=(0, (2, 6)))
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
        ax.minorticks_on()
    axes[-1].set_xlabel("S [mm]")
    axes[0].set_title("MSOL step-0 electric field components at witness positions")
    axes[0].legend(loc="best")
    fig.tight_layout(rect=(0.0, 0.10, 1.0, 0.98))
    add_footer(fig, directory, config_path)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_step0_field_xy_contour(
    trajectory: pd.DataFrame,
    output: Path,
    directory: Path,
    config: dict,
    step: int = 0,
    levels: np.ndarray | None = None,
    grid_data: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    track12 = load_track12_module()

    step_data = trajectory.loc[trajectory["step"] == step].copy()
    if step_data.empty:
        raise ValueError(f"trajectory has no step-{step} rows")
    time_s = float(step_data["time_s"].iloc[0])
    time_ps = 1.0e12 * time_s
    extent_um = config["field_xy_extent_um"]
    points = config["field_xy_points"]
    if grid_data is None:
        x_um, y_um, e_dot_e = evaluate_xy_edot(config, time_s)
    else:
        x_um, y_um, e_dot_e = grid_data

    if levels is None:
        positive = e_dot_e[e_dot_e > 0.0]
        levels = np.geomspace(positive.min(), positive.max(), 42) if positive.size else 42
    sigma_x_um = 1.0e6 * track12.SIGMA_XY_M
    sigma_y_um = 1.0e6 * track12.SIGMA_XY_M
    mid = points // 2
    x_profile = e_dot_e[mid, :]
    y_profile = e_dot_e[:, mid]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 8,
            "savefig.dpi": 220,
        }
    )
    fig = plt.figure(figsize=(8.8, 8.0))
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.15, 4.8, 0.28),
        height_ratios=(1.15, 4.8),
        wspace=0.10,
        hspace=0.10,
    )
    ax_top = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax = fig.add_subplot(grid[1, 1])
    cax = fig.add_subplot(grid[1, 2])

    contour = ax.contourf(
        x_um,
        y_um,
        e_dot_e,
        levels=levels,
        cmap="magma",
        norm=LogNorm(vmin=float(levels[0]), vmax=float(levels[-1])),
    )
    cbar = fig.colorbar(contour, cax=cax)
    cbar.set_label(r"$\mathbf{E}\cdot\mathbf{E}$ [(V/m)$^2$]")

    ax_top.semilogy(x_um, x_profile, color="0.15", lw=1.2)
    ax_top.axvline(-sigma_x_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.axvline(+sigma_x_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.set_xlim(-extent_um, extent_um)
    ax_top.set_ylabel(r"$E\cdot E$")
    ax_top.set_title(r"MSOL transverse field: $\mathbf{E}\cdot\mathbf{E}$")
    ax_top.grid(True, which="both", alpha=0.25)
    ax_top.tick_params(labelbottom=False)

    ax_left.semilogx(y_profile, y_um, color="0.15", lw=1.2)
    ax_left.axhline(-sigma_y_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_left.axhline(+sigma_y_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_left.set_ylim(-extent_um, extent_um)
    ax_left.set_xlabel(r"$E\cdot E$")
    ax_left.set_ylabel(r"y [$\mu$m]")
    ax_left.grid(True, which="both", alpha=0.25)

    witness = step_data.drop_duplicates(["species", "pair"]).copy()
    markers = {"electron": "o", "positron": "s"}
    colors = {"electron": ELECTRON_COLOR, "positron": POSITRON_COLOR}
    labels = {"electron": "e-", "positron": "e+"}
    for species, group in witness.groupby("species", sort=True):
        ax.scatter(
            1.0e6 * group["x"],
            1.0e6 * group["y"],
            s=58 if species == "electron" else 28,
            marker=markers.get(species, "o"),
            facecolors="none" if species == "electron" else colors.get(species, "0.2"),
            edgecolors=colors.get(species, "0.2"),
            linewidths=1.2,
            label=labels.get(species, species),
        )

    ax.axvline(-sigma_x_um, color="white", lw=1.0, ls=(0, (4, 4)), label=r"$\pm\sigma_x,\pm\sigma_y$")
    ax.axvline(+sigma_x_um, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axhline(-sigma_y_um, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axhline(+sigma_y_um, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel("")
    ax.set_xlim(-extent_um, extent_um)
    ax.set_ylim(-extent_um, extent_um)
    ax.text(
        0.02,
        0.98,
        f"step {step}\nt = {time_ps:.4g} ps",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.45", "alpha": 0.86},
    )
    ax.legend(loc="upper right")
    ax.minorticks_on()
    fig.subplots_adjust(bottom=0.12, top=0.95, left=0.10, right=0.94)
    add_footer(fig, directory, config["config_path"])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_step0_field_zr_contour(
    trajectory: pd.DataFrame,
    output: Path,
    directory: Path,
    config: dict,
    step: int = 0,
    levels: np.ndarray | None = None,
    grid_data: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    track12 = load_track12_module()

    step_data = trajectory.loc[trajectory["step"] == step].copy()
    if step_data.empty:
        raise ValueError(f"trajectory has no step-{step} rows")
    time_s = float(step_data["time_s"].iloc[0])
    time_ps = 1.0e12 * time_s
    z_extent_um = config["field_z_extent_um"]
    r_extent_um = config["field_r_extent_um"]
    points = config["field_zr_points"]
    if grid_data is None:
        z_um, r_um, e_dot_e = evaluate_zr_edot(config, time_s)
    else:
        z_um, r_um, e_dot_e = grid_data

    if levels is None:
        positive = e_dot_e[e_dot_e > 0.0]
        levels = np.geomspace(positive.min(), positive.max(), 42) if positive.size else 42
    sigma_z_um = 1.0e6 * config["source_sigma_z_m"]
    sigma_r_um = 1.0e6 * track12.SIGMA_XY_M
    z_mid = points // 2
    r_axis = 0
    z_profile = e_dot_e[r_axis, :]
    r_profile = e_dot_e[:, z_mid]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 8,
            "savefig.dpi": 220,
        }
    )
    fig = plt.figure(figsize=(9.4, 7.6))
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.15, 5.4, 0.28),
        height_ratios=(1.15, 4.6),
        wspace=0.10,
        hspace=0.10,
    )
    ax_top = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax = fig.add_subplot(grid[1, 1])
    cax = fig.add_subplot(grid[1, 2])

    contour = ax.contourf(
        z_um,
        r_um,
        e_dot_e,
        levels=levels,
        cmap="magma",
        norm=LogNorm(vmin=float(levels[0]), vmax=float(levels[-1])),
    )
    cbar = fig.colorbar(contour, cax=cax)
    cbar.set_label(r"$\mathbf{E}\cdot\mathbf{E}$ [(V/m)$^2$]")

    ax_top.semilogy(z_um, z_profile, color="0.15", lw=1.2)
    ax_top.axvline(-sigma_z_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.axvline(+sigma_z_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.set_xlim(-z_extent_um, z_extent_um)
    ax_top.set_ylabel(r"$E\cdot E$")
    ax_top.set_title(r"MSOL longitudinal-radial field: $\mathbf{E}\cdot\mathbf{E}$")
    ax_top.grid(True, which="both", alpha=0.25)
    ax_top.tick_params(labelbottom=False)

    ax_left.semilogx(r_profile, r_um, color="0.15", lw=1.2)
    ax_left.axhline(sigma_r_um, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_left.set_ylim(0.0, r_extent_um)
    ax_left.set_xlabel(r"$E\cdot E$")
    ax_left.set_ylabel(r"r [$\mu$m]")
    ax_left.grid(True, which="both", alpha=0.25)

    witness = step_data.drop_duplicates(["species", "pair"]).copy()
    witness["z_um"] = 1.0e6 * witness["s"]
    witness["r_um"] = 1.0e6 * np.sqrt(witness["x"] * witness["x"] + witness["y"] * witness["y"])
    markers = {"electron": "o", "positron": "s"}
    colors = {"electron": ELECTRON_COLOR, "positron": POSITRON_COLOR}
    labels = {"electron": "e-", "positron": "e+"}
    for species, group in witness.groupby("species", sort=True):
        ax.scatter(
            group["z_um"],
            group["r_um"],
            s=58 if species == "electron" else 28,
            marker=markers.get(species, "o"),
            facecolors="none" if species == "electron" else colors.get(species, "0.2"),
            edgecolors=colors.get(species, "0.2"),
            linewidths=1.2,
            label=labels.get(species, species),
        )

    ax.axvline(-sigma_z_um, color="white", lw=1.0, ls=(0, (4, 4)), label=r"$\pm\sigma_z,\sigma_r$")
    ax.axvline(+sigma_z_um, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axhline(sigma_r_um, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.set_xlabel(r"z [$\mu$m]")
    ax.set_ylabel("")
    ax.set_xlim(-z_extent_um, z_extent_um)
    ax.set_ylim(0.0, r_extent_um)
    ax.text(
        0.02,
        0.98,
        f"step {step}\nt = {time_ps:.4g} ps",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.45", "alpha": 0.86},
    )
    ax.legend(loc="upper right")
    ax.minorticks_on()
    fig.subplots_adjust(bottom=0.12, top=0.95, left=0.10, right=0.94)
    add_footer(fig, directory, config["config_path"])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def write_contour_sequences(trajectory: pd.DataFrame, directory: Path, config: dict) -> list[Path]:
    count = config["field_contour_steps"]
    stride = config["field_contour_stride"]
    if stride <= 0:
        raise ValueError("field_contour_stride must be positive")

    steps = sorted(int(step) for step in trajectory["step"].drop_duplicates())
    selected_steps = steps[::stride]
    if count > 0:
        selected_steps = selected_steps[:count]
    if not selected_steps:
        return []
    outputs: list[Path] = []
    xy_dir = config["field_xy_sequence_dir"]
    zr_dir = config["field_zr_sequence_dir"]

    if xy_dir is not None:
        xy_dir.mkdir(parents=True, exist_ok=True)
        xy_frames = [
            evaluate_xy_edot(config, float(trajectory.loc[trajectory["step"] == step, "time_s"].iloc[0]))
            for step in selected_steps
        ]
        xy_levels = edot_levels([frame[2] for frame in xy_frames])
        for step, frame in zip(selected_steps, xy_frames, strict=True):
            output = xy_dir / f"msol-edot-xy-step-{step:04d}.png"
            plot_step0_field_xy_contour(
                trajectory,
                output,
                directory,
                config,
                step=step,
                levels=xy_levels,
                grid_data=frame,
            )
            outputs.append(output)

    if zr_dir is not None:
        zr_dir.mkdir(parents=True, exist_ok=True)
        zr_frames = [
            evaluate_zr_edot(config, float(trajectory.loc[trajectory["step"] == step, "time_s"].iloc[0]))
            for step in selected_steps
        ]
        zr_levels = edot_levels([frame[2] for frame in zr_frames])
        for step, frame in zip(selected_steps, zr_frames, strict=True):
            output = zr_dir / f"msol-edot-zr-step-{step:04d}.png"
            plot_step0_field_zr_contour(
                trajectory,
                output,
                directory,
                config,
                step=step,
                levels=zr_levels,
                grid_data=frame,
            )
            outputs.append(output)

    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    config["figure_dir"].mkdir(parents=True, exist_ok=True)
    config["data_dir"].mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(config["cain_csv"])
    manufactured = pd.read_csv(config["trajectory_csv"])

    trajectory_plot = config["trajectory_plot"] or config["figure_dir"] / "cain-fit-half-charge-analytic-vs-cain.png"
    field_plot = config["field_plot"] or config["figure_dir"] / "cain-fit-half-charge-field-vs-time.png"
    step0_field_plot = config["field_components_step0_plot"]
    step0_xy_contour_plot = config["field_xy_step0_contour_plot"]
    step0_zr_contour_plot = config["field_zr_step0_contour_plot"]
    field_summary = config["field_summary_csv"] or config["data_dir"] / "cain-fit-half-charge-field-vs-time-summary.csv"
    directory = Path.cwd().resolve()
    plot_trajectory_comparison(
        reference,
        manufactured,
        trajectory_plot,
        directory,
        config["config_path"],
        config["source_selection"],
        config["charge_scale"],
    )
    field_summary = plot_field_vs_time(
        manufactured,
        field_plot,
        field_summary,
        directory,
        config["config_path"],
    )
    if step0_field_plot is not None:
        plot_step0_field_components(manufactured, step0_field_plot, directory, config["config_path"])
    if step0_xy_contour_plot is not None:
        plot_step0_field_xy_contour(manufactured, step0_xy_contour_plot, directory, config)
    if step0_zr_contour_plot is not None:
        plot_step0_field_zr_contour(manufactured, step0_zr_contour_plot, directory, config)
    sequence_outputs = write_contour_sequences(manufactured, directory, config)

    print(f"Wrote {trajectory_plot}")
    print(f"Wrote {field_plot}")
    print(f"Wrote {field_summary}")
    if step0_field_plot is not None:
        print(f"Wrote {step0_field_plot}")
    if step0_xy_contour_plot is not None:
        print(f"Wrote {step0_xy_contour_plot}")
    if step0_zr_contour_plot is not None:
        print(f"Wrote {step0_zr_contour_plot}")
    if sequence_outputs:
        print(f"Wrote {len(sequence_outputs)} contour sequence plots")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
