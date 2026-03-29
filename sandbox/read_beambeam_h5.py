#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def load_h5py():
    try:
        import h5py  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "h5py is not installed. Activate the local environment first:\n"
            "  source .venv-h5/bin/activate"
        ) from exc
    return h5py


def load_numpy():
    try:
        import numpy as np  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit("numpy is not installed.") from exc
    return np


def load_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "matplotlib is not installed. Install it in the local environment first:\n"
            "  source .venv-h5/bin/activate\n"
            "  python -m pip install matplotlib"
        ) from exc
    return plt


def decode_value(value):
    if hasattr(value, "shape"):
        try:
            if value.shape == ():
                value = value[()]
            elif value.shape == (1,):
                value = value[0]
            else:
                return [decode_value(item) for item in value]
        except Exception:
            pass

    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        try:
            return decode_value(value.item())
        except Exception:
            pass
    return value


def format_scalar(value) -> str:
    value = decode_value(value)
    if isinstance(value, float):
        return f"{value:.9g}"
    if isinstance(value, int):
        return str(value)
    return str(value)


def format_vector(value) -> str:
    value = decode_value(value)
    if not isinstance(value, list):
        return format_scalar(value)
    parts = []
    for item in value:
        if isinstance(item, float):
            parts.append(f"{item:.6g}")
        else:
            parts.append(str(item))
    return "(" + ", ".join(parts) + ")"


def translate_snapshot(value) -> str:
    translated = {
        "normal_tracking": "normal tracking",
        "before_beambeam_mesh_enlarge": "normal tracking",
        "active_beambeam": "beambeam tracking",
    }
    return translated.get(value, value)


def collect_step_rows(h5file):
    rows = []
    step_names = sorted(h5file.keys(), key=step_sort_key)

    for step_name in step_names:
        step = h5file[step_name]
        attrs = step.attrs
        path_length_s = decode_value(attrs.get("path_length_s", "?"))
        shape = decode_value(attrs.get("shape", "?"))
        mesh_origin = decode_value(attrs.get("mesh_origin", "?"))
        mesh_spacing = decode_value(attrs.get("mesh_spacing", "?"))

        snapshot_kind = format_scalar(attrs.get("snapshot_kind", "?"))
        rows.append(
            {
                "step": format_scalar(attrs.get("global_step", "?")),
                "state": translate_snapshot(snapshot_kind),
                "state_raw": snapshot_kind,
                "t": format_scalar(attrs.get("time", "?")),
                "s": format_scalar(attrs.get("path_length_s", "?")),
                "<s>": format_scalar(attrs.get("particle_mean_s", "?")),
                "IP": format_scalar(attrs.get("interaction_point_s", "?")),
                "extend (lab)": format_vector(attrs.get("ip_element_s_range", "?")),
                "Q": format_scalar(attrs.get("particle_total_charge", "?")),
                "shape": format_vector(attrs.get("shape", "?")),
                "origin_z (beam)": format_scalar(mesh_origin[2]) if isinstance(mesh_origin, list) and len(mesh_origin) >= 3 else "?",
                "spacing": format_vector(attrs.get("mesh_spacing", "?")),
            }
        )

    return step_names, rows


def select_states(step_names, rows, state_raw):
    if state_raw is None:
        preferred_states = (
            "active_beambeam",
            "normal_tracking",
        )
        selected_step_names = []
        selected_rows = []
        by_step = {}
        order = {}
        for step_name, row in zip(step_names, rows):
            try:
                step = int(row["step"])
            except ValueError:
                continue
            by_step.setdefault(step, []).append((step_name, row))
            order.setdefault(step, len(order))

        for step in sorted(by_step.keys(), key=lambda x: order[x]):
            entries = by_step[step]
            selected = None
            for preferred_state in preferred_states:
                for step_name, row in entries:
                    if row["state_raw"] == preferred_state:
                        selected = (step_name, row)
                        break
                if selected is not None:
                    break
            if selected is None:
                selected = entries[0]
            selected_step_names.append(selected[0])
            selected_rows.append(selected[1])

        return selected_step_names, selected_rows

    filtered_step_names = []
    filtered_rows = []
    for step_name, row in zip(step_names, rows):
        if row["state_raw"] == state_raw:
            filtered_step_names.append(step_name)
            filtered_rows.append(row)
    return filtered_step_names, filtered_rows


def filter_steps(step_names, rows, start, end):
    if start is None and end is None:
        return step_names, rows

    filtered_step_names = []
    filtered_rows = []

    for step_name, row in zip(step_names, rows):
        try:
            step = int(row["step"])
        except ValueError:
            continue

        if start is not None and step < start:
            continue
        if end is not None and step > end:
            continue

        filtered_step_names.append(step_name)
        filtered_rows.append(row)

    return filtered_step_names, filtered_rows


def print_table(rows) -> None:
    columns = [
        "step",
        "state",
        "t",
        "s",
        "<s>",
        "IP",
        "extend (lab)",
        "Q",
        "shape",
        "origin_z (beam)",
        "spacing",
    ]
    widths = {
        col: max(len(col), *(len(str(row[col])) for row in rows)) if rows else len(col)
        for col in columns
    }

    header = "  ".join(f"{col:<{widths[col]}}" for col in columns)
    print(header)
    print("  ".join("-" * widths[col] for col in columns))
    for row in rows:
        print("  ".join(f"{str(row[col]):<{widths[col]}}" for col in columns))


def step_sort_key(name: str) -> tuple[int, str]:
    if name.startswith("Step#"):
        suffix = name.split("#", 1)[1]
        try:
            return (int(suffix), name)
        except ValueError:
            pass
    return (sys.maxsize, name)


def get_step_rho(step, np):
    rho_group = step["Block"]["rho"]
    dataset_names = sorted(rho_group.keys(), key=step_sort_key)
    if not dataset_names:
        raise SystemExit("rho dataset missing in HDF5 step.")
    rho = np.asarray(rho_group[dataset_names[0]][()])
    origin = np.asarray(rho_group.attrs["__Origin__"], dtype=float)
    spacing = np.asarray(rho_group.attrs["__Spacing__"], dtype=float)
    return rho, origin, spacing


def compute_line_density_data(step, row, bins, np):
    rho, origin, spacing = get_step_rho(step, np)

    # H5hut exposes the 3D block as (z, y, x) for this writer layout.
    nz = rho.shape[0]
    z_centers = origin[2] + (np.arange(nz) + 0.5) * spacing[2]

    line_density_native = rho.sum(axis=(1, 2)) * spacing[0] * spacing[1]

    z_edges = origin[2] + np.arange(nz + 1) * spacing[2]
    native_charge = line_density_native * spacing[2]

    target_edges = np.linspace(z_edges[0], z_edges[-1], bins + 1)
    charge_per_bin, _ = np.histogram(
        z_centers,
        bins=target_edges,
        weights=native_charge,
    )
    target_widths = np.diff(target_edges)
    line_density = charge_per_bin / target_widths
    target_centers = 0.5 * (target_edges[:-1] + target_edges[1:])

    return {
        "step": row["step"],
        "state": row["state"],
        "z": target_centers,
        "line_density": line_density,
    }


def compute_xz_projection_data(step, row, np, side):
    rho, origin, spacing = get_step_rho(step, np)

    nx = rho.shape[2]
    nz = rho.shape[0]
    x0 = origin[0]
    x1 = origin[0] + nx * spacing[0]
    z0 = origin[2]
    z1 = origin[2] + nz * spacing[2]

    attrs = step.attrs
    path_length_s = decode_value(attrs.get("path_length_s", "?"))
    interaction_point_s = decode_value(attrs.get("interaction_point_s", "?"))
    particle_mean_r = decode_value(attrs.get("particle_mean_r", "?"))
    particle_mean_s = decode_value(attrs.get("particle_mean_s", "?"))

    z_shift = path_length_s if isinstance(path_length_s, float) else 0.0

    primary_centroid = None
    secondary_centroid = None
    primary_mask = None
    is_beambeam_tracking = row["state"] == "beambeam tracking"

    if (
        isinstance(particle_mean_r, list)
        and len(particle_mean_r) >= 3
        and isinstance(path_length_s, float)
        and isinstance(interaction_point_s, float)
    ):
        primary_x = float(particle_mean_r[0])
        primary_z_lab = float(path_length_s + particle_mean_r[2])
        primary_centroid = (primary_x, primary_z_lab)
        if is_beambeam_tracking:
            secondary_centroid = (
                primary_x,
                float(2.0 * interaction_point_s - primary_z_lab),
            )
    elif (
        isinstance(path_length_s, float)
        and isinstance(interaction_point_s, float)
        and isinstance(particle_mean_s, float)
    ):
        ip_beam_z = interaction_point_s - path_length_s
        x_centers = origin[0] + (np.arange(nx) + 0.5) * spacing[0]
        z_centers = origin[2] + (np.arange(nz) + 0.5) * spacing[2]

        # Use only the physical negative charge contribution for centroiding.
        weights = np.maximum(-rho, 0.0)

        if particle_mean_s <= interaction_point_s:
            z_mask = z_centers[:, None, None] <= ip_beam_z
        else:
            z_mask = z_centers[:, None, None] >= ip_beam_z

        primary_mask = z_mask
        primary_weights = weights * z_mask
        primary_weight_sum = float(primary_weights.sum())
        if primary_weight_sum > 0.0:
            x_grid = x_centers[None, None, :]
            z_grid = z_centers[:, None, None]
            primary_x = float((primary_weights * x_grid).sum() / primary_weight_sum)
            primary_z_beam = float((primary_weights * z_grid).sum() / primary_weight_sum)
            primary_centroid = (primary_x, primary_z_beam)
            if is_beambeam_tracking:
                secondary_centroid = (primary_x, 2.0 * ip_beam_z - primary_z_beam)

    if side == "combined" or primary_mask is None:
        rho_for_projection = rho
    elif side == "primary":
        rho_for_projection = np.where(primary_mask, rho, 0.0)
    elif side == "secondary":
        rho_for_projection = np.where(primary_mask, 0.0, rho)
    else:
        raise SystemExit(f"unsupported gallery side: {side}")

    # H5hut exposes the 3D block as (z, y, x) for this writer layout.
    projection_xz = rho_for_projection.sum(axis=1) * spacing[1]

    if primary_centroid is not None:
        if not (
            isinstance(particle_mean_r, list)
            and len(particle_mean_r) >= 3
            and isinstance(path_length_s, float)
            and isinstance(interaction_point_s, float)
        ):
            primary_centroid = (primary_centroid[0], primary_centroid[1] + z_shift)
    if secondary_centroid is not None:
        if not (
            isinstance(particle_mean_r, list)
            and len(particle_mean_r) >= 3
            and isinstance(path_length_s, float)
            and isinstance(interaction_point_s, float)
        ):
            secondary_centroid = (secondary_centroid[0], secondary_centroid[1] + z_shift)

    return {
        "step": row["step"],
        "state": row["state"],
        "projection": projection_xz,
        "extent": [x0, x1, z0 + z_shift, z1 + z_shift],
        "primary_centroid": primary_centroid,
        "secondary_centroid": secondary_centroid,
    }


def plot_line_density_z(path: Path, bins: int, start, end, gif_path, state_raw) -> int:
    if bins <= 0:
        raise SystemExit("--bins must be > 0")

    h5py = load_h5py()
    np = load_numpy()
    plt = load_matplotlib()
    animation = None
    if gif_path is not None:
        from matplotlib import animation as mpl_animation  # type: ignore
        animation = mpl_animation

    with h5py.File(path, "r") as h5file:
        step_names, rows = collect_step_rows(h5file)
        step_names, rows = select_states(step_names, rows, state_raw)
        step_names, rows = filter_steps(step_names, rows, start, end)
        if not step_names:
            raise SystemExit("no steps found in selected range")

        line_density_steps = []
        y_min = None
        y_max = None
        x_min = None
        x_max = None

        for step_name, row in zip(step_names, rows):
            step = h5file[step_name]
            data = compute_line_density_data(step, row, bins, np)
            line_density_steps.append(data)

            current_y_min = float(data["line_density"].min())
            current_y_max = float(data["line_density"].max())
            current_x_min = float(data["z"][0])
            current_x_max = float(data["z"][-1])

            y_min = current_y_min if y_min is None else min(y_min, current_y_min)
            y_max = current_y_max if y_max is None else max(y_max, current_y_max)
            x_min = current_x_min if x_min is None else min(x_min, current_x_min)
            x_max = current_x_max if x_max is None else max(x_max, current_x_max)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_xlabel("z [m]")
        ax.set_ylabel("line density [C/m]")
        ax.set_title(f"BeamBeam line density in z ({bins} bins)")
        ax.grid(True, alpha=0.25)
        ax.set_xlim(x_min, x_max)
        if y_min is not None and y_max is not None:
            if y_min == y_max:
                pad = 1.0 if y_min == 0.0 else 0.05 * abs(y_min)
            else:
                pad = 0.05 * (y_max - y_min)
            ax.set_ylim(y_min - pad, y_max + pad)

        if gif_path is None:
            for data in line_density_steps:
                label = f"step {data['step']}"
                ax.step(
                    data["z"],
                    data["line_density"],
                    where="mid",
                    linewidth=1.0,
                    alpha=0.8,
                    label=label,
                )
            ax.legend(fontsize=8, ncol=2)
            fig.tight_layout()
            plt.show()
            return 0

        line, = ax.step([], [], where="mid", linewidth=2.0)
        title = ax.set_title("")

        def update(frame_index):
            data = line_density_steps[frame_index]
            line.set_data(data["z"], data["line_density"])
            title.set_text(
                f"BeamBeam line density in z ({bins} bins) | step {data['step']} | {data['state']}"
            )
            return line, title

        anim = animation.FuncAnimation(
            fig,
            update,
            frames=len(line_density_steps),
            interval=300,
            blit=False,
            repeat=True,
        )
        fig.tight_layout()
        anim.save(str(gif_path), writer="pillow", fps=4)
        plt.close(fig)

    return 0


def plot_gallery_xz(path: Path, start, end, cols: int, side: str, state_raw) -> int:
    if cols <= 0:
        raise SystemExit("--gallery-cols must be > 0")

    h5py = load_h5py()
    np = load_numpy()
    plt = load_matplotlib()

    with h5py.File(path, "r") as h5file:
        step_names, rows = collect_step_rows(h5file)
        step_names, rows = select_states(step_names, rows, state_raw)
        step_names, rows = filter_steps(step_names, rows, start, end)
        if not step_names:
            raise SystemExit("no steps found in selected range")

        panels = []
        vmin = None
        vmax = None

        for step_name, row in zip(step_names, rows):
            step = h5file[step_name]
            data = compute_xz_projection_data(step, row, np, side)
            panels.append(data)

            current_min = float(data["projection"].min())
            current_max = float(data["projection"].max())
            vmin = current_min if vmin is None else min(vmin, current_min)
            vmax = current_max if vmax is None else max(vmax, current_max)

        if vmin is None or vmax is None:
            raise SystemExit("failed to determine gallery color scale")

        if vmin < 0.0 and vmax <= 0.0:
            color_limits = (vmin, 0.0)
        elif vmin >= 0.0 and vmax > 0.0:
            color_limits = (0.0, vmax)
        else:
            bound = max(abs(vmin), abs(vmax))
            color_limits = (-bound, bound)

        n_panels = len(panels)
        n_cols = min(cols, n_panels)
        n_rows = (n_panels + n_cols - 1) // n_cols

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(4.0 * n_cols, 3.2 * n_rows),
            squeeze=False,
        )
        fig.patch.set_facecolor("white")

        for ax, panel in zip(axes.flat, panels):
            x_limits = (panel["extent"][0], panel["extent"][1])
            z_limits = (panel["extent"][2], panel["extent"][3])
            ax.imshow(
                panel["projection"],
                origin="lower",
                aspect="auto",
                extent=panel["extent"],
                vmin=color_limits[0],
                vmax=color_limits[1],
                cmap="bone_r",
            )
            if panel["primary_centroid"] is not None:
                ax.plot(
                    panel["primary_centroid"][0],
                    panel["primary_centroid"][1],
                    marker="x",
                    color="black",
                    markersize=8,
                    markeredgewidth=1.5,
                )
            if panel["secondary_centroid"] is not None:
                ax.plot(
                    panel["secondary_centroid"][0],
                    panel["secondary_centroid"][1],
                    marker="o",
                    markerfacecolor="none",
                    markeredgecolor="black",
                    markersize=8,
                    markeredgewidth=1.5,
                )
            ax.set_xlim(*x_limits)
            ax.set_ylim(*z_limits)
            ax.set_title(f"step {panel['step']}", fontsize=10, pad=4)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)

        for ax in axes.flat[n_panels:]:
            ax.axis("off")

        side_label = {
            "combined": "combined",
            "primary": "primary side only",
            "secondary": "mirrored side only",
        }[side]
        fig.suptitle(
            f"BeamBeam projected charge density in x,z ({side_label})",
            fontsize=12,
        )
        fig.tight_layout()
        plt.show()

    return 0


def overview(path: Path, table: bool, start, end) -> int:
    h5py = load_h5py()

    with h5py.File(path, "r") as h5file:
        step_names, rows = collect_step_rows(h5file)
        step_names, rows = filter_steps(step_names, rows, start, end)

        print(f"file: {path}")
        print(f"steps: {len(step_names)}")
        print()

        if table:
            print_table(rows)
            return 0

        for step_name, row in zip(step_names, rows):
            print(f"{step_name}:")
            print(f"  step: {row['step']}")
            print(f"  state: {row['state']}")
            print(f"  t: {row['t']}")
            print(f"  s: {row['s']}")
            print(f"  <s>: {row['<s>']}")
            print(f"  IP: {row['IP']}")
            print(f"  extend (lab): {row['extend (lab)']}")
            print(f"  Q: {row['Q']}")
            print(f"  shape: {row['shape']}")
            print(f"  origin_z (beam): {row['origin_z (beam)']}")
            print(f"  spacing: {row['spacing']}")
            print()

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read OPALX BeamBeam HDF5 diagnostics."
    )
    parser.add_argument("filename", help="BeamBeam diagnostics HDF5 file")
    parser.add_argument(
        "--overview",
        action="store_true",
        help="Show step metadata overview for all steps",
    )
    parser.add_argument(
        "--table",
        action="store_true",
        help="Render overview output as one table row per step",
    )
    parser.add_argument(
        "--line-density-z",
        action="store_true",
        help="Plot line density in z for every stored step",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=32,
        help="Number of z bins for --line-density-z",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=None,
        help="First global step to include",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="Last global step to include",
    )
    parser.add_argument(
        "--gif",
        default=None,
        help="Write the line-density-z plot as an animated GIF",
    )
    parser.add_argument(
        "--gallery-xz",
        action="store_true",
        help="Show an x,z projected charge-density gallery for all selected steps",
    )
    parser.add_argument(
        "--gallery-cols",
        type=int,
        default=4,
        help="Number of columns for --gallery-xz",
    )
    parser.add_argument(
        "--gallery-xz-side",
        choices=("combined", "primary", "secondary"),
        default="combined",
        help="Select which side of the BeamBeam density to show in --gallery-xz",
    )
    parser.add_argument(
        "--state",
        choices=(
            "normal_tracking",
            "active_beambeam",
            "before_beambeam_mesh_enlarge",
        ),
        default=None,
        help="Select one stored BeamBeam diagnostic state. Default for plots is one entry per step, preferring active_beambeam.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    path = Path(args.filename)
    if not path.exists():
        raise SystemExit(f"file not found: {path}")
    if args.start is not None and args.end is not None and args.start > args.end:
        raise SystemExit("--start must be <= --end")

    if args.overview:
        return overview(path, args.table, args.start, args.end)
    if args.line_density_z:
        gif_path = Path(args.gif) if args.gif is not None else None
        return plot_line_density_z(
            path,
            args.bins,
            args.start,
            args.end,
            gif_path,
            args.state,
        )
    if args.gallery_xz:
        return plot_gallery_xz(
            path,
            args.start,
            args.end,
            args.gallery_cols,
            args.gallery_xz_side,
            args.state,
        )

    parser.error("no action selected; use --overview, --line-density-z, or --gallery-xz")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
