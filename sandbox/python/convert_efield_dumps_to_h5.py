#!/usr/bin/env python3
"""Convert BeamBeam EF_vector debug dumps to HDF5 and render an x-z gallery."""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = ROOT / "sandbox/track-e-p/data"
DEFAULT_OUTPUT_H5 = ROOT / "sandbox/data/gamma_gamma_large_cylinder_efield_debug.h5"
DEFAULT_GALLERY = ROOT / "sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug.png"
DEFAULT_EZ_GALLERY = ROOT / "sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug_ez.png"
DEFAULT_TIMELINE = ROOT / "sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug_lab_timeline.png"
DEFAULT_PROFILE = ROOT / "sandbox/note/figs/gamma_gamma_large_cylinder_efield_debug_z_profile.png"
HEADER_RE = re.compile(r"([A-Za-z0-9_.\[\]]+)\s*=\s*(\([^)]*\)|\[[^\]]*\]|[^#\s]+)")
STEP_RE = re.compile(r"-(\d{6})\.dat$")


def parse_numeric(text: str) -> float:
    return float(text.strip().strip("[]"))


def parse_vector(text: str) -> tuple[float, ...]:
    values = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", text)
    return tuple(float(value) for value in values)


def step_from_path(path: Path) -> int:
    match = STEP_RE.search(path.name)
    if match is None:
        return -1
    return int(match.group(1))


def parse_ef_dump(path: Path) -> dict[str, object]:
    metadata: dict[str, str] = {}
    rows: list[tuple[int, int, int, float, float, float, float, float, float]] = []

    with path.open() as stream:
        for line in stream:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                for key, value in HEADER_RE.findall(line[1:].strip()):
                    metadata[key] = value.strip()
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            rows.append(
                (
                    int(parts[0]),
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    float(parts[4]),
                    float(parts[5]),
                    float(parts[6]),
                    float(parts[7]),
                    float(parts[8]),
                )
            )

    if not rows:
        raise ValueError(f"no EF samples found in {path}")

    data = np.asarray(rows, dtype=float)
    i_vals = data[:, 0].astype(int)
    j_vals = data[:, 1].astype(int)
    k_vals = data[:, 2].astype(int)
    imin, jmin, kmin = i_vals.min(), j_vals.min(), k_vals.min()
    nx = i_vals.max() - imin + 1
    ny = j_vals.max() - jmin + 1
    nz = k_vals.max() - kmin + 1

    x = np.full(nx, np.nan)
    y = np.full(ny, np.nan)
    z = np.full(nz, np.nan)
    ex = np.zeros((nx, ny, nz), dtype=np.float64)
    ey = np.zeros_like(ex)
    ez = np.zeros_like(ex)

    for i_raw, j_raw, k_raw, x_val, y_val, z_val, ex_val, ey_val, ez_val in rows:
        i = i_raw - imin
        j = j_raw - jmin
        k = k_raw - kmin
        x[i] = x_val
        y[j] = y_val
        z[k] = z_val
        ex[i, j, k] = ex_val
        ey[i, j, k] = ey_val
        ez[i, j, k] = ez_val

    if np.isnan(x).any() or np.isnan(y).any() or np.isnan(z).any():
        raise ValueError(f"incomplete EF grid in {path}")

    return {
        "path": path,
        "metadata": metadata,
        "x": x,
        "y": y,
        "z": z,
        "Ex": ex,
        "Ey": ey,
        "Ez": ez,
    }


def metadata_float(metadata: dict[str, str], key: str) -> float:
    value = metadata.get(key)
    if value is None:
        return math.nan
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    return parse_numeric(value)


def metadata_vector(metadata: dict[str, str], key: str) -> tuple[float, ...]:
    value = metadata.get(key)
    if value is None:
        return ()
    return parse_vector(value)


def select_paths(args: argparse.Namespace) -> list[Path]:
    pattern = f"{args.stem}-EF_vector-beambeam_e-*.dat"
    paths = sorted(args.input_dir.glob(pattern), key=step_from_path)
    if args.start_step is not None:
        paths = [path for path in paths if step_from_path(path) >= args.start_step]
    if args.end_step is not None:
        paths = [path for path in paths if step_from_path(path) <= args.end_step]
    if args.stride > 1:
        paths = paths[:: args.stride]
    if args.max_frames is not None and len(paths) > args.max_frames:
        indices = np.linspace(0, len(paths) - 1, args.max_frames).round().astype(int)
        paths = [paths[int(index)] for index in indices]
    if not paths:
        raise FileNotFoundError(f"no EF_vector dumps found for {pattern} in {args.input_dir}")
    return paths


def write_h5(paths: list[Path], output: Path) -> list[dict[str, object]]:
    output.parent.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    with h5py.File(output, "w") as h5:
        h5.attrs["source"] = "OPALX BeamBeam EF_vector debug dumps"
        for ordinal, path in enumerate(paths):
            dump = parse_ef_dump(path)
            metadata = dump["metadata"]
            group = h5.create_group(f"Step#{step_from_path(path):06d}")
            group.attrs["source_file"] = str(path)
            for key, value in metadata.items():
                group.attrs[key] = value
            for axis in ("x", "y", "z"):
                group.create_dataset(axis, data=dump[axis], compression="gzip")
            block = group.create_group("E")
            for component in ("Ex", "Ey", "Ez"):
                block.create_dataset(component, data=dump[component], compression="gzip")
            e_abs = np.sqrt(dump["Ex"] ** 2 + dump["Ey"] ** 2 + dump["Ez"] ** 2)
            block.create_dataset("Eabs", data=e_abs, compression="gzip")

            mean_r = metadata_vector(metadata, "particle_mean_r")
            ip_local_z = metadata_float(metadata, "interaction_point_local_z")
            path_length_s = metadata_float(metadata, "path_length_s")
            ip_s = metadata_float(metadata, "interaction_point_s")
            primary_z = mean_r[2] if len(mean_r) >= 3 else math.nan
            mirror_z = 2.0 * ip_local_z - primary_z if math.isfinite(ip_local_z) else math.nan
            primary_lab_z = path_length_s + primary_z if math.isfinite(path_length_s) else math.nan
            mirror_lab_z = path_length_s + mirror_z if math.isfinite(path_length_s) else math.nan
            group.attrs["primary_local_z_m"] = primary_z
            group.attrs["mirror_local_z_m"] = mirror_z
            group.attrs["primary_lab_z_m"] = primary_lab_z
            group.attrs["mirror_lab_z_m"] = mirror_lab_z
            group.attrs["interaction_point_lab_z_m"] = ip_s
            summaries.append(
                {
                    "ordinal": ordinal,
                    "step": step_from_path(path),
                    "global_step": metadata_float(metadata, "global_step"),
                    "time_ps": 1.0e12 * metadata_float(metadata, "time"),
                    "ip_local_z": ip_local_z,
                    "ip_lab_z": ip_s,
                    "primary_z": primary_z,
                    "mirror_z": mirror_z,
                    "primary_lab_z": primary_lab_z,
                    "mirror_lab_z": mirror_lab_z,
                    "path": path,
                }
            )
    return summaries


def render_gallery(h5_path: Path, output: Path, max_frames: int, component: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(output.parent / ".plot-cache/matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output.parent / ".plot-cache/xdg"))
    import matplotlib.pyplot as plt

    with h5py.File(h5_path, "r") as h5:
        step_names = sorted(h5.keys(), key=lambda name: int(name.split("#", 1)[1]))
        if len(step_names) > max_frames:
            indices = np.linspace(0, len(step_names) - 1, max_frames).round().astype(int)
            step_names = [step_names[int(index)] for index in indices]

        cols = min(4, len(step_names))
        rows = math.ceil(len(step_names) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.25 * rows), squeeze=False)
        values: list[np.ndarray] = []
        payloads = []
        for name in step_names:
            group = h5[name]
            x = np.asarray(group["x"])
            z = np.asarray(group["z"])
            field = np.asarray(group[f"E/{component}"])
            mid_y = field.shape[1] // 2
            plane = field[:, mid_y, :].T
            values.append(plane)
            payloads.append((name, group, x, z, plane))

        if component == "Eabs":
            vmin = 0.0
            vmax = max(float(np.nanpercentile(value, 99.0)) for value in values)
            vmax = vmax if vmax > 0 else 1.0
            cmap = "magma"
            colorbar_label = r"$|E|$"
        else:
            vmax = max(float(np.nanpercentile(np.abs(value), 99.0)) for value in values)
            vmax = vmax if vmax > 0 else 1.0
            vmin = -vmax
            cmap = "coolwarm"
            colorbar_label = rf"${component}$"

        last_image = None
        for ax, item in zip(axes.flat, payloads, strict=False):
            name, group, x, z, plane = item
            extent = [float(z[0]), float(z[-1]), float(x[0]), float(x[-1])]
            last_image = ax.imshow(
                plane,
                origin="lower",
                extent=extent,
                aspect="auto",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )
            ip_local_z = metadata_float(dict(group.attrs), "interaction_point_local_z")
            mean_r = metadata_vector(dict(group.attrs), "particle_mean_r")
            primary_z = mean_r[2] if len(mean_r) >= 3 else math.nan
            mirror_z = 2.0 * ip_local_z - primary_z if math.isfinite(ip_local_z) else math.nan
            if math.isfinite(ip_local_z):
                ax.axvline(ip_local_z, color="white", linewidth=1.1, label="local IP")
            if math.isfinite(primary_z):
                ax.axvline(primary_z, color="#2AA6B8", linewidth=1.0, linestyle="-", label="bunch[0] local")
            if math.isfinite(mirror_z):
                ax.axvline(mirror_z, color="#D33682", linewidth=1.0, linestyle="--", label="mirror local")
            time_ps = metadata_float(dict(group.attrs), "time") * 1.0e12
            ax.set_title(f"{name}, t={time_ps:.1f} ps")
            ax.set_xlabel("local z [m]")
            ax.set_ylabel("x [m]")
            ax.grid(color="white", alpha=0.18, linewidth=0.5)
            if ax is axes.flat[0]:
                ax.legend(loc="upper right", fontsize=7, framealpha=0.75)
        for ax in axes.flat[len(payloads) :]:
            ax.axis("off")
        fig.subplots_adjust(left=0.06, right=0.88, bottom=0.06, top=0.91, hspace=0.45, wspace=0.30)
        if last_image is not None:
            cbar_ax = fig.add_axes((0.91, 0.12, 0.018, 0.74))
            fig.colorbar(last_image, cax=cbar_ax, label=colorbar_label)
        fig.suptitle(
            f"BeamBeam mirrored E-field debug: {component}, central y slice in local coordinates",
            y=0.975,
        )
        fig.savefig(output, dpi=180)
        plt.close(fig)


def render_lab_timeline(summaries: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(output.parent / ".plot-cache/matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output.parent / ".plot-cache/xdg"))
    import matplotlib.pyplot as plt

    time_ps = np.asarray([row["time_ps"] for row in summaries], dtype=float)
    primary_lab_z = np.asarray([row["primary_lab_z"] for row in summaries], dtype=float)
    mirror_lab_z = np.asarray([row["mirror_lab_z"] for row in summaries], dtype=float)
    ip_lab_z = np.asarray([row["ip_lab_z"] for row in summaries], dtype=float)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(time_ps, primary_lab_z, color="#2AA6B8", linewidth=2.0, label="bunch[0] lab z")
    ax.plot(time_ps, mirror_lab_z, color="#D33682", linewidth=2.0, linestyle="--", label="mirrored lab z")
    ax.plot(time_ps, ip_lab_z, color="0.2", linewidth=1.4, label="IP lab z")
    ax.set_xlabel("time [ps]")
    ax.set_ylabel("global z / s [m]")
    ax.set_title("BeamBeam source and mirrored source positions in the lab frame")
    ax.grid(alpha=0.28, linewidth=0.7)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def render_z_profile_gallery(h5_path: Path, output: Path, max_frames: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(output.parent / ".plot-cache/matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output.parent / ".plot-cache/xdg"))
    import matplotlib.pyplot as plt

    with h5py.File(h5_path, "r") as h5:
        step_names = sorted(h5.keys(), key=lambda name: int(name.split("#", 1)[1]))
        if len(step_names) > max_frames:
            indices = np.linspace(0, len(step_names) - 1, max_frames).round().astype(int)
            step_names = [step_names[int(index)] for index in indices]

        cols = min(4, len(step_names))
        rows = math.ceil(len(step_names) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 2.8 * rows), squeeze=False)
        for ax, name in zip(axes.flat, step_names, strict=False):
            group = h5[name]
            z = np.asarray(group["z"])
            e_abs = np.asarray(group["E/Eabs"])
            profile = e_abs.sum(axis=(0, 1))
            if np.nanmax(profile) > 0.0:
                profile = profile / np.nanmax(profile)

            ip_local_z = metadata_float(dict(group.attrs), "interaction_point_local_z")
            primary_z = metadata_float(dict(group.attrs), "primary_local_z_m")
            mirror_z = metadata_float(dict(group.attrs), "mirror_local_z_m")
            time_ps = metadata_float(dict(group.attrs), "time") * 1.0e12

            ax.plot(z, profile, color="0.12", linewidth=1.8)
            if math.isfinite(primary_z):
                ax.axvline(primary_z, color="#2AA6B8", linewidth=1.0, label="bunch[0] local")
            if math.isfinite(ip_local_z):
                ax.axvline(ip_local_z, color="0.45", linewidth=1.0, label="local IP")
            if math.isfinite(mirror_z):
                ax.axvline(mirror_z, color="#D33682", linewidth=1.0, linestyle="--", label="mirror local")
            ax.set_title(f"{name}, t={time_ps:.1f} ps")
            ax.set_xlabel("local z [m]")
            ax.set_ylabel(r"$\sum_{x,y}|E|$ / max")
            ax.set_ylim(0.0, 1.08)
            ax.grid(alpha=0.28, linewidth=0.6)
            if ax is axes.flat[0]:
                ax.legend(loc="upper right", fontsize=7, framealpha=0.8)
        for ax in axes.flat[len(step_names) :]:
            ax.axis("off")
        fig.suptitle("BeamBeam E-field longitudinal profile, transverse integral", y=0.98)
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))
        fig.savefig(output, dpi=180)
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--stem", default="gamma_gamma_pairs-large-cylinder")
    parser.add_argument("--output-h5", type=Path, default=DEFAULT_OUTPUT_H5)
    parser.add_argument("--gallery", type=Path, default=DEFAULT_GALLERY)
    parser.add_argument("--ez-gallery", type=Path, default=DEFAULT_EZ_GALLERY)
    parser.add_argument("--timeline", type=Path, default=DEFAULT_TIMELINE)
    parser.add_argument("--z-profile", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--start-step", type=int)
    parser.add_argument("--end-step", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int)
    parser.add_argument("--gallery-frames", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paths = select_paths(args)
    summaries = write_h5(paths, args.output_h5)
    render_gallery(args.output_h5, args.gallery, args.gallery_frames, "Eabs")
    render_gallery(args.output_h5, args.ez_gallery, args.gallery_frames, "Ez")
    render_lab_timeline(summaries, args.timeline)
    render_z_profile_gallery(args.output_h5, args.z_profile, args.gallery_frames)
    print(f"read {len(paths)} EF dumps")
    print(f"wrote {args.output_h5}")
    print(f"wrote {args.gallery}")
    print(f"wrote {args.ez_gallery}")
    print(f"wrote {args.timeline}")
    print(f"wrote {args.z_profile}")
    for row in (summaries[0], summaries[len(summaries) // 2], summaries[-1]):
        print(
            "step={step} time={time_ps:.3f} ps local: ip_z={ip_local_z:.6g} "
            "primary_z={primary_z:.6g} mirror_z={mirror_z:.6g}; lab: "
            "ip_z={ip_lab_z:.6g} primary_z={primary_lab_z:.6g} "
            "mirror_z={mirror_lab_z:.6g}".format(**row)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
