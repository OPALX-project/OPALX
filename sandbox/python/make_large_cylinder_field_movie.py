#!/usr/bin/env python3
"""Animate the large-cylinder BeamBeam crossing histogram with source fields.

The particle positions are read from the OPALX H5 trajectory.  The electric
field background is a vectorized version of the spherical boosted Gaussian
source-field model used by ``sandbox/python/boosted_gaussian_witness.py``.  This
is a diagnostic overlay, not an OPALX field dump: it shows the expected primary
source Coulomb field while the source is active, then switches the overlay off
after the BeamBeam source retirement time.

The movie is written to ``sandbox/data`` by default, while the preview PNG is
written to ``sandbox/note/figs`` so the note can use it without depending on the
old regression directory.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from pathlib import Path

import h5py
import imageio.v2 as imageio
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TRACK_DIR = ROOT / "sandbox/track-e-p"
DEFAULT_OUTPUT = ROOT / "sandbox/data/gamma_gamma_large_cylinder_field_histogram.mp4"
DEFAULT_PREVIEW = ROOT / "sandbox/note/figs/gamma_gamma_large_cylinder_field_histogram_preview.png"

C_LIGHT = 299_792_458.0
EPSILON_0 = 8.854_187_812_8e-12
ELEMENTARY_CHARGE_C = 1.602_176_634e-19
ELECTRON_REST_EV = 510_998.95


@dataclass(frozen=True)
class SpeciesData:
    label: str
    color: str
    times_s: np.ndarray
    x_by_step: list[np.ndarray]
    y_by_step: list[np.ndarray]
    z_by_step: list[np.ndarray]
    first_crossing_z: np.ndarray
    first_crossing_t: np.ndarray


def configure_environment(output: Path) -> None:
    cache_dir = ROOT / "sandbox/data/.plot-cache"
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))
    output.parent.mkdir(parents=True, exist_ok=True)


def beta_from_kinetic_energy_mev(kinetic_mev: float) -> float:
    gamma = 1.0 + kinetic_mev * 1.0e6 / ELECTRON_REST_EV
    return math.sqrt(1.0 - 1.0 / (gamma * gamma))


def load_species(stem: str, container: str, label: str, color: str, aperture_radius_m: float) -> SpeciesData:
    path = TRACK_DIR / f"{stem}_{container}.h5"
    times: list[float] = []
    x_by_step: list[np.ndarray] = []
    y_by_step: list[np.ndarray] = []
    z_by_step: list[np.ndarray] = []
    first_z: dict[int, float] = {}
    first_t: dict[int, float] = {}

    with h5py.File(path, "r") as h5:
        steps = sorted(
                (key for key in h5.keys() if key.startswith("Step#")),
                key=lambda key: int(key.split("#", 1)[1]))
        for step in steps:
            group = h5[step]
            time_s = float(group.attrs["TIME"][0])
            ids = np.asarray(group["id"], dtype=int)
            x = np.asarray(group["x"], dtype=float)
            y = np.asarray(group["y"], dtype=float)
            z = np.asarray(group["z"], dtype=float)
            radius = np.hypot(x, y)
            outside = radius >= aperture_radius_m
            for particle_id, z_position in zip(ids[outside], z[outside], strict=True):
                first_z.setdefault(int(particle_id), float(z_position))
                first_t.setdefault(int(particle_id), time_s)
            times.append(time_s)
            x_by_step.append(x)
            y_by_step.append(y)
            z_by_step.append(z)

    ordered_ids = sorted(first_z)
    return SpeciesData(
            label=label,
            color=color,
            times_s=np.asarray(times, dtype=float),
            x_by_step=x_by_step,
            y_by_step=y_by_step,
            z_by_step=z_by_step,
            first_crossing_z=np.asarray([first_z[particle_id] for particle_id in ordered_ids], dtype=float),
            first_crossing_t=np.asarray([first_t[particle_id] for particle_id in ordered_ids], dtype=float),
    )


def active_frame_indices(n_steps: int, max_frames: int) -> np.ndarray:
    early = np.arange(0, min(n_steps, 42), 2, dtype=int)
    late = np.arange(min(n_steps, 50), n_steps, 12, dtype=int)
    indices = np.unique(np.concatenate([early, late, np.asarray([n_steps - 1], dtype=int)]))
    if len(indices) <= max_frames:
        return indices
    selected = np.linspace(0, len(indices) - 1, max_frames).round().astype(int)
    return np.unique(indices[selected])


def gaussian_source_field_xz(
        z_grid: np.ndarray,
        x_grid: np.ndarray,
        time_s: float,
        source_charge_c: float,
        source_sigma_m: float,
        source_beta: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Ex, Ez, and |E| for a boosted spherical Gaussian source."""
    gamma = 1.0 / math.sqrt(1.0 - source_beta * source_beta)
    source_z = source_beta * C_LIGHT * time_s
    dx = x_grid
    dz = z_grid - source_z
    r_prime_x = dx
    r_prime_z = gamma * dz
    radius2 = r_prime_x * r_prime_x + r_prime_z * r_prime_z
    radius = np.sqrt(radius2)

    coeff = source_charge_c / (4.0 * math.pi * EPSILON_0)
    sigma2 = source_sigma_m * source_sigma_m
    safe_radius = np.maximum(radius, 1.0e-18)
    u = safe_radius / (math.sqrt(2.0) * source_sigma_m)
    erf_u = np.vectorize(math.erf)(u)
    exp_term = np.exp(-0.5 * radius2 / sigma2)
    shape = erf_u - math.sqrt(2.0 / math.pi) * (safe_radius / source_sigma_m) * exp_term
    scale = coeff * shape / (safe_radius**3)
    scale = np.where(radius > 0.0, scale, 0.0)

    e_prime_x = scale * r_prime_x
    e_prime_z = scale * r_prime_z
    e_lab_x = gamma * e_prime_x
    e_lab_z = e_prime_z
    e_abs = np.hypot(e_lab_x, e_lab_z)
    return e_lab_x, e_lab_z, e_abs


def make_movie(args: argparse.Namespace) -> None:
    configure_environment(args.output)
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 160,
        "font.size": 10,
    })

    electron = load_species(args.stem, "c1", r"$e^-$", "#2AA6B8", args.aperture_radius_m)
    positron = load_species(args.stem, "c2", r"$e^+$", "#D33682", args.aperture_radius_m)
    display_label = args.label if args.label is not None else args.stem
    if len(electron.times_s) != len(positron.times_s):
        raise ValueError("electron and positron H5 files have different step counts")
    times_s = electron.times_s
    frame_indices = active_frame_indices(len(times_s), args.max_frames)

    z_axis = np.linspace(args.z_min_m, args.z_max_m, args.field_grid_z)
    x_axis = np.linspace(-args.aperture_radius_m, args.aperture_radius_m, args.field_grid_x)
    z_grid, x_grid = np.meshgrid(z_axis, x_axis)
    qz_axis = np.linspace(args.cylinder_begin_m, args.cylinder_end_m, args.quiver_z)
    qx_axis = np.linspace(-0.85 * args.aperture_radius_m, 0.85 * args.aperture_radius_m, args.quiver_x)
    qz_grid, qx_grid = np.meshgrid(qz_axis, qx_axis)

    source_beta = beta_from_kinetic_energy_mev(args.primary_kinetic_mev)
    source_charge_c = -args.primary_electrons * ELEMENTARY_CHARGE_C

    active_logs = []
    for time_s in times_s[frame_indices]:
        if time_s <= args.retire_time_s:
            _ex, _ez, e_abs = gaussian_source_field_xz(
                    z_grid, x_grid, time_s, source_charge_c, args.source_sigma_m, source_beta)
            active_logs.append(np.log10(np.maximum(e_abs, args.field_floor_v_per_m)))
    if active_logs:
        stacked = np.concatenate([values.ravel() for values in active_logs])
        vmin = float(np.percentile(stacked, 3.0))
        vmax = float(np.percentile(stacked, 99.0))
    else:
        vmin, vmax = -4.0, 4.0
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin >= vmax:
        vmin, vmax = -4.0, 4.0

    edges = np.linspace(args.hist_min_m, args.hist_max_m, args.hist_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist_width = 0.42 * (edges[1] - edges[0])

    fig, (ax_field, ax_hist) = plt.subplots(
            2, 1, figsize=(10.0, 7.0), height_ratios=(1.45, 1.0), constrained_layout=True)

    initial_log = np.full_like(z_grid, vmin)
    image = ax_field.imshow(
            initial_log,
            origin="lower",
            extent=(args.z_min_m, args.z_max_m, -args.aperture_radius_m, args.aperture_radius_m),
            aspect="auto",
            cmap="magma",
            vmin=vmin,
            vmax=vmax)
    colorbar = fig.colorbar(image, ax=ax_field, pad=0.012)
    colorbar.set_label(r"$\log_{10}(|E| / [\mathrm{V/m}])$")
    quiver = ax_field.quiver(qz_grid, qx_grid, np.zeros_like(qz_grid), np.zeros_like(qx_grid),
                             color="white", pivot="mid", scale=34.0, width=0.0025, alpha=0.86)
    electron_points = ax_field.scatter([], [], s=6, color=electron.color, alpha=0.60, edgecolors="none", label=r"$e^-$")
    positron_points = ax_field.scatter([], [], s=6, color=positron.color, alpha=0.60, edgecolors="none", label=r"$e^+$")
    source_marker = ax_field.scatter([], [], s=70, marker="x", color="white", linewidths=1.6, label="primary center")
    title = ax_field.set_title("")
    ax_field.axhline(args.aperture_radius_m, color="white", linewidth=0.8, alpha=0.75)
    ax_field.axhline(-args.aperture_radius_m, color="white", linewidth=0.8, alpha=0.75)
    ax_field.axvline(args.cylinder_begin_m, color="white", linewidth=0.8, linestyle=":", alpha=0.75)
    ax_field.axvline(args.cylinder_end_m, color="white", linewidth=0.8, linestyle=":", alpha=0.75)
    ax_field.axvline(0.5 * (args.cylinder_begin_m + args.cylinder_end_m), color="white",
                     linewidth=0.8, linestyle="--", alpha=0.75)
    ax_field.text(
            0.5 * (args.cylinder_begin_m + args.cylinder_end_m),
            -0.92 * args.aperture_radius_m,
            "IP",
            color="white",
            ha="center",
            va="top")
    ax_field.set_xlim(args.z_min_m, args.z_max_m)
    ax_field.set_ylim(-args.aperture_radius_m, args.aperture_radius_m)
    ax_field.set_ylabel("x [m]")
    ax_field.legend(loc="upper right", frameon=True, fontsize=8)

    e_bars = ax_hist.bar(centers - 0.5 * hist_width, np.zeros_like(centers), width=hist_width,
                         color=electron.color, edgecolor="black", linewidth=0.35, alpha=0.78, label=r"$e^-$")
    p_bars = ax_hist.bar(centers + 0.5 * hist_width, np.zeros_like(centers), width=hist_width,
                         color=positron.color, edgecolor="black", linewidth=0.35, alpha=0.78, label=r"$e^+$")
    ax_hist.axvspan(args.cylinder_begin_m, args.cylinder_end_m, color="0.92", zorder=-2)
    ax_hist.axvline(0.5 * (args.cylinder_begin_m + args.cylinder_end_m), color="0.2", linewidth=0.8)
    ax_hist.set_xlim(args.hist_min_m, args.hist_max_m)
    ax_hist.set_ylim(0, args.hist_ymax)
    ax_hist.set_xlabel("first cylinder-edge crossing z [m]")
    ax_hist.set_ylabel("cumulative first crossings")
    ax_hist.grid(axis="y", color="0.86", linewidth=0.6)
    ax_hist.spines[["top", "right"]].set_visible(False)
    ax_hist.legend(loc="upper right", frameon=False, ncols=2)
    hist_title = ax_hist.set_title("")

    def particle_offsets(species: SpeciesData, step_index: int) -> np.ndarray:
        z = species.z_by_step[step_index]
        x = species.x_by_step[step_index]
        return np.column_stack([z, x])

    def update_frame(step_index: int) -> None:
        time_s = float(times_s[step_index])
        time_ps = 1.0e12 * time_s
        source_active = time_s <= args.retire_time_s

        if source_active:
            ex, ez, e_abs = gaussian_source_field_xz(
                    z_grid, x_grid, time_s, source_charge_c, args.source_sigma_m, source_beta)
            log_e = np.log10(np.maximum(e_abs, args.field_floor_v_per_m))
            qex, qez, qe_abs = gaussian_source_field_xz(
                    qz_grid, qx_grid, time_s, source_charge_c, args.source_sigma_m, source_beta)
            qnorm = np.maximum(qe_abs, args.field_floor_v_per_m)
            quiver.set_UVC(qez / qnorm, qex / qnorm)
            source_z = source_beta * C_LIGHT * time_s
            source_marker.set_offsets([[source_z, 0.0]])
            field_state = "analytic primary field active"
        else:
            log_e = np.full_like(z_grid, vmin)
            quiver.set_UVC(np.zeros_like(qz_grid), np.zeros_like(qx_grid))
            source_marker.set_offsets(np.empty((0, 2)))
            field_state = "source retired: field overlay set to zero"

        image.set_data(log_e)
        electron_points.set_offsets(particle_offsets(electron, step_index))
        positron_points.set_offsets(particle_offsets(positron, step_index))

        e_counts, _ = np.histogram(
                electron.first_crossing_z[electron.first_crossing_t <= time_s], bins=edges)
        p_counts, _ = np.histogram(
                positron.first_crossing_z[positron.first_crossing_t <= time_s], bins=edges)
        for patch, height in zip(e_bars, e_counts, strict=True):
            patch.set_height(float(height))
        for patch, height in zip(p_bars, p_counts, strict=True):
            patch.set_height(float(height))

        title.set_text(
                f"{display_label}: t = {time_ps:.1f} ps, "
                f"retire = {1.0e12 * args.retire_time_s:.1f} ps | {field_state}")
        hist_title.set_text(
                f"crossings so far: e- {int(e_counts.sum())}/{len(electron.first_crossing_z)}, "
                f"e+ {int(p_counts.sum())}/{len(positron.first_crossing_z)}")

    writer_kwargs = {"fps": args.fps}
    if args.output.suffix.lower() == ".mp4":
        writer_kwargs.update({"codec": "libx264", "quality": 8, "pixelformat": "yuv420p"})

    with imageio.get_writer(args.output, **writer_kwargs) as writer:
        preview_written = False
        for step_index in frame_indices:
            update_frame(int(step_index))
            fig.canvas.draw()
            frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
            writer.append_data(frame)
            if not preview_written and times_s[step_index] >= args.preview_time_ps * 1.0e-12:
                fig.savefig(args.preview, dpi=160)
                preview_written = True
        if not preview_written:
            fig.savefig(args.preview, dpi=160)
    plt.close(fig)

    print(f"wrote {args.output}")
    print(f"wrote {args.preview}")
    print(f"frames: {len(frame_indices)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", default="gamma_gamma_pairs-large-cylinder")
    parser.add_argument("--label", default=None, help="Short display label for the movie title.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--preview", type=Path, default=DEFAULT_PREVIEW)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--max-frames", type=int, default=110)
    parser.add_argument("--preview-time-ps", type=float, default=1120.0)
    parser.add_argument("--aperture-radius-m", type=float, default=0.15)
    parser.add_argument("--cylinder-begin-m", type=float, default=0.01)
    parser.add_argument("--cylinder-end-m", type=float, default=0.33)
    parser.add_argument("--z-min-m", type=float, default=-0.06)
    parser.add_argument("--z-max-m", type=float, default=0.47)
    parser.add_argument("--hist-min-m", type=float, default=-0.06)
    parser.add_argument("--hist-max-m", type=float, default=0.47)
    parser.add_argument("--hist-bins", type=int, default=18)
    parser.add_argument("--hist-ymax", type=float, default=80.0)
    parser.add_argument("--field-grid-z", type=int, default=150)
    parser.add_argument("--field-grid-x", type=int, default=86)
    parser.add_argument("--quiver-z", type=int, default=19)
    parser.add_argument("--quiver-x", type=int, default=9)
    parser.add_argument("--field-floor-v-per-m", type=float, default=1.0e-6)
    parser.add_argument("--primary-kinetic-mev", type=float, default=245.0)
    parser.add_argument("--primary-electrons", type=float, default=1.25e10)
    parser.add_argument("--source-sigma-m", type=float, default=0.6e-3)
    parser.add_argument("--retire-time-s", type=float, default=571.35e-12)
    return parser.parse_args()


def main() -> int:
    make_movie(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
