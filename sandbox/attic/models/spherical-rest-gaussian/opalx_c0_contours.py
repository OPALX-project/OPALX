#!/usr/bin/env python3
"""Make OPALX C0 particle-sampled E dot E contour plots."""

from __future__ import annotations

import argparse
import json
import os
import textwrap
from pathlib import Path

import h5py
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPT_DIR / "opalx_c0_contours.json"


def resolve(config_dir: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return config_dir / path


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config_dir = path.resolve().parent
    config["config_path"] = path.resolve()
    config["h5_file"] = resolve(config_dir, config["h5_file"])
    config["figure_dir"] = resolve(config_dir, config["figure_dir"])
    config["xy_sequence_dir"] = resolve(config_dir, config["xy_sequence_dir"])
    config["zr_sequence_dir"] = resolve(config_dir, config["zr_sequence_dir"])
    return config


def configure_matplotlib() -> None:
    cache_dir = SCRIPT_DIR / "atic" / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)


def add_footer(fig, directory: Path, config_path: Path, h5_file: Path) -> None:
    wrap_width = max(80, int(fig.get_figwidth() * 14))
    footer_text = "\n".join(
        [
            textwrap.fill(f"directory: {directory}", width=wrap_width),
            textwrap.fill(f"config: {config_path}", width=wrap_width),
            textwrap.fill(f"h5: {h5_file}", width=wrap_width),
        ]
    )
    fig.text(0.5, 0.01, footer_text, ha="center", va="bottom", fontsize=5.8, color="0.35")


def sorted_step_keys(handle: h5py.File) -> list[str]:
    return sorted(handle.keys(), key=lambda key: int(key.split("#", maxsplit=1)[1]))


def bin_mean_2d(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> np.ndarray:
    sums, _, _ = np.histogram2d(y, x, bins=(y_edges, x_edges), weights=values)
    counts, _, _ = np.histogram2d(y, x, bins=(y_edges, x_edges))
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = sums / counts
    mean[counts == 0.0] = np.nan
    return mean


def bin_center(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def load_step_grid(group: h5py.Group, config: dict) -> dict:
    x_um = np.asarray(group["x"], dtype=float) * 1.0e6
    y_um = np.asarray(group["y"], dtype=float) * 1.0e6
    z_um_abs = np.asarray(group["z"], dtype=float) * 1.0e6
    z0_um = float(np.mean(z_um_abs))
    z_source_um = z_um_abs - z0_um
    z_reference_um = float(config.get("zr_z_reference_um", z0_um))
    z_um = z_um_abs - z_reference_um
    r_um = np.sqrt(x_um * x_um + y_um * y_um)
    ex = np.asarray(group["Ex"], dtype=float)
    ey = np.asarray(group["Ey"], dtype=float)
    ez = np.asarray(group["Ez"], dtype=float)
    edot = ex * ex + ey * ey + ez * ez
    bbox_incr = float(config.get("bbox_increment_percent", 1.0)) / 100.0
    mesh_bounds = {}
    for axis_name, values in (("x", x_um), ("y", y_um), ("z", z_um_abs)):
        axis_min = float(np.min(values))
        axis_max = float(np.max(values))
        axis_span = axis_max - axis_min
        mesh_bounds[axis_name] = (
            axis_min - bbox_incr * axis_span,
            axis_max + bbox_incr * axis_span,
        )

    xy_extent = float(config["xy_extent_um"])
    xy_edges = np.linspace(-xy_extent, xy_extent, int(config["xy_bins"]) + 1)
    z_half_width = float(config["xy_z_slice_half_width_um"])
    xy_mask = np.abs(z_source_um) <= z_half_width
    if not np.any(xy_mask):
        raise ValueError(f"empty x-y z-slice for {group.name}")
    xy_grid = bin_mean_2d(x_um[xy_mask], y_um[xy_mask], edot[xy_mask], xy_edges, xy_edges)

    z_extent = float(config["zr_z_extent_um"])
    r_extent = float(config["zr_r_extent_um"])
    z_edges = np.linspace(-z_extent, z_extent, int(config["zr_bins_z"]) + 1)
    r_edges = np.linspace(0.0, r_extent, int(config["zr_bins_r"]) + 1)
    zr_mask = (np.abs(z_um) <= z_extent) & (r_um <= r_extent)
    if not np.any(zr_mask):
        raise ValueError(f"empty z-r range for {group.name}")
    zr_grid = bin_mean_2d(z_um[zr_mask], r_um[zr_mask], edot[zr_mask], z_edges, r_edges)

    return {
        "step": int(group.attrs["Step"][0]),
        "global_step": int(group.attrs["GlobalTrackStep"][0]),
        "time_s": float(group.attrs["TIME"][0]),
        "total_count": int(len(x_um)),
        "z0_um": z0_um,
        "z_reference_um": z_reference_um,
        "mesh_bounds_um": mesh_bounds,
        "xy_edges": xy_edges,
        "xy_grid": xy_grid,
        "xy_count": int(np.count_nonzero(xy_mask)),
        "z_edges": z_edges,
        "r_edges": r_edges,
        "zr_grid": zr_grid,
        "zr_count": int(np.count_nonzero(zr_mask)),
    }


def load_witness_frame(config: dict, step: int) -> dict[str, np.ndarray]:
    h5_file = Path(config["h5_file"])
    stem = h5_file.name
    if "_c0" not in stem:
        return {}
    species_files = {
        "electron": h5_file.with_name(stem.replace("_c0", "_c1")),
        "positron": h5_file.with_name(stem.replace("_c0", "_c2")),
    }
    z_reference_um = float(config.get("zr_z_reference_um", 0.0))
    output: dict[str, np.ndarray] = {}
    group_name = f"Step#{step}"
    for species, path in species_files.items():
        if not path.exists():
            continue
        with h5py.File(path, "r") as handle:
            if group_name not in handle:
                continue
            group = handle[group_name]
            x_um = np.asarray(group["x"], dtype=float) * 1.0e6
            y_um = np.asarray(group["y"], dtype=float) * 1.0e6
            z_um = np.asarray(group["z"], dtype=float) * 1.0e6 - z_reference_um
            r_um = np.sqrt(x_um * x_um + y_um * y_um)
            ids = np.asarray(group["id"], dtype=int)
            output[species] = np.column_stack((ids + 1, z_um, r_um))
    return output


def positive_limits(grids: list[np.ndarray]) -> tuple[float, float]:
    positives = [grid[np.isfinite(grid) & (grid > 0.0)] for grid in grids]
    positive = np.concatenate([part for part in positives if part.size])
    if not positive.size:
        raise ValueError("no positive E dot E values found")
    return float(np.nanmin(positive)), float(np.nanmax(positive))


def nanmean_axis(values: np.ndarray, axis: int) -> np.ndarray:
    valid = np.isfinite(values)
    sums = np.nansum(values, axis=axis)
    counts = np.sum(valid, axis=axis)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = sums / counts
    mean[counts == 0] = np.nan
    return mean


def plot_xy(frame: dict, config: dict, output: Path, vmin: float, vmax: float) -> None:
    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "savefig.dpi": 220,
        }
    )
    x_edges = frame["xy_edges"]
    y_edges = frame["xy_edges"]
    x_centers = bin_center(x_edges)
    y_centers = bin_center(y_edges)
    grid = np.ma.masked_invalid(frame["xy_grid"])
    x_profile = nanmean_axis(frame["xy_grid"], axis=0)
    y_profile = nanmean_axis(frame["xy_grid"], axis=1)
    extent = float(config["xy_extent_um"])

    fig = plt.figure(figsize=(8.8, 8.0))
    layout = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.15, 4.8, 0.28),
        height_ratios=(1.15, 4.8),
        wspace=0.10,
        hspace=0.10,
    )
    ax_top = fig.add_subplot(layout[0, 1])
    ax_left = fig.add_subplot(layout[1, 0])
    ax = fig.add_subplot(layout[1, 1])
    cax = fig.add_subplot(layout[1, 2])

    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        grid,
        cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
    )
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(r"$\mathbf{E}\cdot\mathbf{E}$ [(V/m)$^2$]")
    cbar.ax.tick_params(labelsize=9)

    sigma_x = float(config["sigma_x_um"])
    sigma_y = float(config["sigma_y_um"])
    ax_top.semilogy(x_centers, x_profile, color="0.15", lw=1.2)
    ax_top.axvline(-sigma_x, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.axvline(+sigma_x, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.set_xlim(-extent, extent)
    ax_top.set_ylabel(r"$E\cdot E$")
    ax_top.set_title(r"OPALX C0 transverse field: $\mathbf{E}\cdot\mathbf{E}$")
    ax_top.grid(True, which="both", alpha=0.25)
    ax_top.tick_params(labelbottom=False)

    ax_left.semilogx(y_profile, y_centers, color="0.15", lw=1.2)
    ax_left.axhline(-sigma_y, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_left.axhline(+sigma_y, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_left.set_ylim(-extent, extent)
    ax_left.set_xlabel(r"$E\cdot E$")
    ax_left.set_ylabel(r"y [$\mu$m]")
    ax_left.grid(True, which="both", alpha=0.25)

    ax.axvline(-sigma_x, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axvline(+sigma_x, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axhline(-sigma_y, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axhline(+sigma_y, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel("")
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.text(
        0.02,
        0.98,
        f"step {frame['step']}\nglobal {frame['global_step']}\nt = {frame['time_s'] * 1.0e12:.4g} ps\n"
        f"|z-<z>| < {config['xy_z_slice_half_width_um']} um\n"
        f"N shown/total = {frame['xy_count']}/{frame['total_count']}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.45", "alpha": 0.86},
    )
    fig.subplots_adjust(bottom=0.14, top=0.95, left=0.10, right=0.90)
    add_footer(fig, Path.cwd().resolve(), config["config_path"], config["h5_file"])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_zr(frame: dict, config: dict, output: Path, vmin: float, vmax: float) -> None:
    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 11,
            "axes.labelsize": 14,
            "axes.titlesize": 14,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "savefig.dpi": 220,
        }
    )
    z_edges = frame["z_edges"]
    r_edges = frame["r_edges"]
    z_centers = bin_center(z_edges)
    r_centers = bin_center(r_edges)
    grid = np.ma.masked_invalid(frame["zr_grid"])
    z_profile = nanmean_axis(frame["zr_grid"], axis=0)
    r_profile = nanmean_axis(frame["zr_grid"], axis=1)
    z_extent = float(config["zr_z_extent_um"])
    r_extent = float(config["zr_r_extent_um"])

    fig = plt.figure(figsize=(9.4, 7.6))
    layout = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.15, 5.4, 0.28),
        height_ratios=(1.15, 4.6),
        wspace=0.10,
        hspace=0.10,
    )
    ax_top = fig.add_subplot(layout[0, 1])
    ax_left = fig.add_subplot(layout[1, 0])
    ax = fig.add_subplot(layout[1, 1])
    cax = fig.add_subplot(layout[1, 2])

    mesh = ax.pcolormesh(
        z_edges,
        r_edges,
        grid,
        cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
    )
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(r"$\mathbf{E}\cdot\mathbf{E}$ [(V/m)$^2$]")
    cbar.ax.tick_params(labelsize=9)

    sigma_z = float(config["sigma_z_um"])
    sigma_r = float(config["sigma_x_um"])
    mesh_bounds = frame["mesh_bounds_um"]
    mesh_z_min = mesh_bounds["z"][0] - frame["z_reference_um"]
    mesh_z_max = mesh_bounds["z"][1] - frame["z_reference_um"]
    mesh_r_extent = min(
        max(abs(mesh_bounds["x"][0]), abs(mesh_bounds["x"][1])),
        max(abs(mesh_bounds["y"][0]), abs(mesh_bounds["y"][1])),
    )
    ax_top.semilogy(z_centers, z_profile, color="0.15", lw=1.2)
    ax_top.axvline(-sigma_z, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.axvline(+sigma_z, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_top.axvline(mesh_z_min, color="#00a6c8", lw=1.0, ls=(0, (6, 3)))
    ax_top.axvline(mesh_z_max, color="#00a6c8", lw=1.0, ls=(0, (6, 3)))
    ax_top.set_xlim(-z_extent, z_extent)
    ax_top.set_ylabel(r"$E\cdot E$")
    ax_top.set_title(r"OPALX C0 longitudinal-radial field: $\mathbf{E}\cdot\mathbf{E}$")
    ax_top.grid(True, which="both", alpha=0.25)
    ax_top.tick_params(labelbottom=False)

    ax_left.semilogx(r_profile, r_centers, color="0.15", lw=1.2)
    ax_left.axhline(sigma_r, color="0.45", lw=1.0, ls=(0, (4, 4)))
    ax_left.axhline(mesh_r_extent, color="#00a6c8", lw=1.0, ls=(0, (6, 3)))
    ax_left.set_ylim(0.0, r_extent)
    ax_left.set_xlabel(r"$E\cdot E$")
    ax_left.set_ylabel(r"r [$\mu$m]")
    ax_left.grid(True, which="both", alpha=0.25)

    ax.axvline(-sigma_z, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axvline(+sigma_z, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axhline(sigma_r, color="white", lw=1.0, ls=(0, (4, 4)))
    ax.axvspan(mesh_z_min, mesh_z_max, ymin=0.0, ymax=min(mesh_r_extent / r_extent, 1.0),
               facecolor="none", edgecolor="#00a6c8", linewidth=1.2, linestyle=(0, (6, 3)),
               label="OPALX mesh extent")
    ax.axhline(mesh_r_extent, color="#00a6c8", lw=1.2, ls=(0, (6, 3)))
    witnesses = load_witness_frame(config, frame["step"])
    witness_styles = {
        "electron": {"marker": "o", "edgecolors": "#c00000", "facecolors": "none", "label": "e- witness"},
        "positron": {"marker": "s", "edgecolors": "#0017b8", "facecolors": "none", "label": "e+ witness"},
    }
    for species, points in witnesses.items():
        style = witness_styles[species]
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
        for pair, z_value, r_value in points:
            dy = 0.12 if species == "electron" else -0.18
            ax.annotate(
                str(int(pair)),
                (z_value, r_value),
                xytext=(3, dy * 10),
                textcoords="offset points",
                fontsize=7,
                color=style["edgecolors"],
                zorder=6,
            )
    z_label = str(config.get("zr_z_label", "z - <z>"))
    ax.set_xlabel(rf"{z_label} [$\mu$m]")
    ax.set_ylabel("")
    ax.set_xlim(-z_extent, z_extent)
    ax.set_ylim(0.0, r_extent)
    ax.text(
        0.02,
        0.98,
        f"step {frame['step']}\nglobal {frame['global_step']}\nt = {frame['time_s'] * 1.0e12:.4g} ps\n"
        f"<z> = {frame['z0_um']:.2f} um\nz_ref = {frame['z_reference_um']:.2f} um\n"
        f"N shown/total = {frame['zr_count']}/{frame['total_count']}\n"
        f"mesh z-IP = [{mesh_z_min:.1f}, {mesh_z_max:.1f}] um",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.45", "alpha": 0.86},
    )
    ax.legend(loc="upper right", fontsize=7, framealpha=0.86)
    fig.subplots_adjust(bottom=0.14, top=0.95, left=0.10, right=0.90)
    add_footer(fig, Path.cwd().resolve(), config["config_path"], config["h5_file"])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    max_steps = int(config.get("max_steps", 0))
    with h5py.File(config["h5_file"], "r") as handle:
        keys = sorted_step_keys(handle)
        if max_steps > 0:
            keys = keys[:max_steps]
        frames = [load_step_grid(handle[key], config) for key in keys]

    xy_vmin, xy_vmax = positive_limits([frame["xy_grid"] for frame in frames])
    zr_vmin, zr_vmax = positive_limits([frame["zr_grid"] for frame in frames])

    outputs: list[Path] = []
    for frame in frames:
        xy_output = config["xy_sequence_dir"] / f"opalx-c0-edot-xy-step-{frame['step']:04d}.png"
        zr_output = config["zr_sequence_dir"] / f"opalx-c0-edot-zr-step-{frame['step']:04d}.png"
        plot_xy(frame, config, xy_output, xy_vmin, xy_vmax)
        plot_zr(frame, config, zr_output, zr_vmin, zr_vmax)
        outputs.extend([xy_output, zr_output])

    print(f"Read {config['h5_file']}")
    print(f"Wrote {len(outputs)} OPALX C0 contour plots")
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
