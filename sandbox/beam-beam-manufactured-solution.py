
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

EPSILON_0 = 8.8541878128e-12


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
            "matplotlib is not installed. Activate the local environment first:\n"
            "  source .venv-h5/bin/activate"
        ) from exc
    return plt


def load_pandas():
    try:
        import pandas as pd  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit("pandas is not installed.") from exc
    return pd


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


def step_sort_key(name: str) -> tuple[int, str]:
    if name.startswith("Step#"):
        suffix = name.split("#", 1)[1]
        try:
            return (int(suffix), name)
        except ValueError:
            pass
    return (sys.maxsize, name)


def erf_array(np, values):
    return np.vectorize(math.erf, otypes=[float])(values)


def gaussian_single_fields(np, x, y, z, center, sigma, charge):
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]

    r2 = dx * dx + dy * dy + dz * dz
    r = np.sqrt(r2)
    sigma2 = sigma * sigma
    u = r / (math.sqrt(2.0) * sigma)
    exp_term = np.exp(-0.5 * r2 / sigma2)

    rho = (
        charge
        / (((2.0 * math.pi) ** 1.5) * sigma**3)
        * exp_term
    )

    coeff = charge / (4.0 * math.pi * EPSILON_0)
    erf_u = erf_array(np, u)

    phi = coeff * np.divide(
        erf_u,
        r,
        out=np.full_like(r, math.sqrt(2.0 / math.pi) / sigma),
        where=r > 0.0,
    )

    shape = erf_u - math.sqrt(2.0 / math.pi) * (r / sigma) * exp_term
    prefactor = np.divide(
        shape,
        r * r * r,
        out=np.zeros_like(r),
        where=r > 0.0,
    )

    ex = coeff * prefactor * dx
    ey = coeff * prefactor * dy
    ez = coeff * prefactor * dz

    return rho, phi, ex, ey, ez


def gaussian_pair_fields(np, x, y, z, sigma, charge, centers):
    shape = (z.shape[0], y.shape[1], x.shape[2])
    rho_total = np.zeros(shape, dtype=float)
    phi_total = np.zeros(shape, dtype=float)
    ex_total = np.zeros(shape, dtype=float)
    ey_total = np.zeros(shape, dtype=float)
    ez_total = np.zeros(shape, dtype=float)

    for center in centers:
        rho, phi, ex, ey, ez = gaussian_single_fields(np, x, y, z, center, sigma, charge)
        rho_total += rho
        phi_total += phi
        ex_total += ex
        ey_total += ey
        ez_total += ez

    return {
        "rho": rho_total,
        "phi": phi_total,
        "Ex": ex_total,
        "Ey": ey_total,
        "Ez": ez_total,
    }


def build_grid(np, origin, spacing, shape, cell_centered=True):
    nx = int(shape[0])
    ny = int(shape[1])
    nz = int(shape[2])

    offset = 0.5 if cell_centered else 0.0
    x = origin[0] + (np.arange(nx) + offset) * spacing[0]
    y = origin[1] + (np.arange(ny) + offset) * spacing[1]
    z = origin[2] + (np.arange(nz) + offset) * spacing[2]

    X = x[None, None, :]
    Y = y[None, :, None]
    Z = z[:, None, None]
    return X, Y, Z


def default_grid(args):
    x_span = args.x_span if args.x_span is not None else 6.0 * args.sigma
    y_span = args.y_span if args.y_span is not None else 6.0 * args.sigma
    z_span_default = 2.0 * (args.half_separation + 3.0 * args.sigma)
    z_span = args.z_span if args.z_span is not None else z_span_default

    origin = (
        -0.5 * x_span,
        -0.5 * y_span,
        -0.5 * z_span,
    )
    spacing = (
        x_span / args.nx,
        y_span / args.ny,
        z_span / args.nz,
    )
    shape = (args.nx, args.ny, args.nz)
    return origin, spacing, shape


def select_h5_step(h5file, step_number, state_raw):
    candidates = []
    for step_name in sorted(h5file.keys(), key=step_sort_key):
        step = h5file[step_name]
        attrs = step.attrs
        global_step = decode_value(attrs.get("global_step"))
        if int(global_step) != step_number:
            continue
        snapshot_kind = str(decode_value(attrs.get("snapshot_kind", "")))
        candidates.append((step_name, snapshot_kind))

    if not candidates:
        raise SystemExit(f"no HDF5 step found for global_step={step_number}")

    if state_raw is not None:
        for candidate in candidates:
            if candidate[1] == state_raw:
                return candidate
        available = ", ".join(snapshot for _, snapshot in candidates)
        raise SystemExit(
            f"step {step_number} does not contain state {state_raw!r}; available: {available}"
        )

    preferred = ("active_beambeam", "normal_tracking", "before_beambeam_mesh_enlarge")
    for snapshot_kind in preferred:
        for candidate in candidates:
            if candidate[1] == snapshot_kind:
                return candidate

    return candidates[0]


def read_h5_scalar_field(np, step, name):
    field_group = step["Block"][name]
    dataset_names = sorted(field_group.keys(), key=step_sort_key)
    if not dataset_names:
        raise SystemExit(f"{name} dataset missing in HDF5 step.")
    field = np.asarray(field_group[dataset_names[0]][()])
    origin = np.asarray(field_group.attrs["__Origin__"], dtype=float)
    spacing = np.asarray(field_group.attrs["__Spacing__"], dtype=float)
    return field, origin, spacing


def parse_header_value(text):
    text = text.strip()
    if text.startswith("(") and text.endswith(")"):
        return [float(item.strip()) for item in text[1:-1].split(",") if item.strip()]
    try:
        return float(text)
    except ValueError:
        return text


def read_ascii_field_dump(np, path: Path):
    pd = load_pandas()
    metadata = {}
    data_start = None
    column_count = None

    with path.open("r", encoding="utf-8") as stream:
        for line_number, raw_line in enumerate(stream):
            line = raw_line.rstrip("\n")
            if not line.startswith("#"):
                data_start = line_number
                column_count = len(line.split())
                break

            header = line[1:].strip()
            if header.startswith("origin="):
                match = re.match(r"origin=\s*(\([^)]+\))\s+h=\s*(\([^)]+\))(?:\s+nghosts=(\S+))?", header)
                if match:
                    metadata["mesh_origin"] = parse_header_value(match.group(1))
                    metadata["mesh_spacing"] = parse_header_value(match.group(2))
                    if match.group(3) is not None:
                        metadata["nghosts"] = int(float(match.group(3)))
                continue
            if "=" in header:
                key, value = header.split("=", 1)
                metadata[key.strip()] = parse_header_value(value.strip())

    if data_start is None or column_count is None:
        raise SystemExit(f"no field rows found in {path}")

    if column_count == 7:
        columns = ["i", "j", "k", "x", "y", "z", "value"]
    elif column_count == 9:
        columns = ["i", "j", "k", "x", "y", "z", "x_value", "y_value", "z_value"]
    else:
        raise SystemExit(f"unsupported field row width {column_count} in {path}")

    frame = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=columns,
        engine="python",
    )

    i_values = sorted(frame["i"].unique())
    j_values = sorted(frame["j"].unique())
    k_values = sorted(frame["k"].unique())
    shape = (len(i_values), len(j_values), len(k_values))
    i_map = {value: index for index, value in enumerate(i_values)}
    j_map = {value: index for index, value in enumerate(j_values)}
    k_map = {value: index for index, value in enumerate(k_values)}

    if column_count == 7:
        field = np.zeros((shape[2], shape[1], shape[0]), dtype=float)
        for row in frame.itertuples(index=False):
            field[k_map[row.k], j_map[row.j], i_map[row.i]] = row.value
    else:
        field = {
            "Ex": np.zeros((shape[2], shape[1], shape[0]), dtype=float),
            "Ey": np.zeros((shape[2], shape[1], shape[0]), dtype=float),
            "Ez": np.zeros((shape[2], shape[1], shape[0]), dtype=float),
        }
        for row in frame.itertuples(index=False):
            iz = k_map[row.k]
            iy = j_map[row.j]
            ix = i_map[row.i]
            field["Ex"][iz, iy, ix] = row.x_value
            field["Ey"][iz, iy, ix] = row.y_value
            field["Ez"][iz, iy, ix] = row.z_value

    origin = np.asarray(metadata["mesh_origin"], dtype=float)
    spacing = np.asarray(metadata["mesh_spacing"], dtype=float)
    return {"metadata": metadata, "field": field, "origin": origin, "spacing": spacing, "shape": shape}


def manufactured_setup_from_ascii(np, rho_path: Path, phi_path: Path, e_path: Path):
    rho_dump = read_ascii_field_dump(np, rho_path)
    phi_dump = read_ascii_field_dump(np, phi_path)
    e_dump = read_ascii_field_dump(np, e_path)

    metadata = rho_dump["metadata"]
    mean_r = metadata.get("particle_mean_r")
    if not (isinstance(mean_r, list) and len(mean_r) >= 3):
        raise SystemExit("particle_mean_r missing or invalid in ASCII dump.")

    primary_center = (float(mean_r[0]), float(mean_r[1]), float(mean_r[2]))
    centers = [primary_center]

    snapshot_kind = str(metadata.get("snapshot_kind", ""))
    if snapshot_kind == "active_beambeam_field_diagnostics":
        if "interaction_point_local_z" in metadata:
            ip_local_z = float(metadata["interaction_point_local_z"])
        else:
            ip_local_z = float(metadata["interaction_point_s"]) - float(metadata["path_length_s"])
        centers.append((primary_center[0], primary_center[1], 2.0 * ip_local_z - primary_center[2]))

    return {
        "snapshot_kind": snapshot_kind,
        "origin": rho_dump["origin"],
        "spacing": rho_dump["spacing"],
        "shape": rho_dump["shape"],
        "charge": float(metadata["particle_total_charge"]),
        "interaction_point_local_z": float(metadata.get("interaction_point_local_z", "nan")),
        "centers": centers,
        "opalx": {
            "rho": rho_dump["field"],
            "phi": phi_dump["field"],
            "Ex": e_dump["field"]["Ex"],
            "Ey": e_dump["field"]["Ey"],
            "Ez": e_dump["field"]["Ez"],
        },
    }


def manufactured_setup_from_h5(np, h5_path: Path, step_number: int, state_raw: str | None):
    h5py = load_h5py()
    with h5py.File(h5_path, "r") as h5file:
        step_name, snapshot_kind = select_h5_step(h5file, step_number, state_raw)
        step = h5file[step_name]
        attrs = step.attrs

        rho, origin, spacing = read_h5_scalar_field(np, step, "rho")
        phi, _, _ = read_h5_scalar_field(np, step, "phi")
        ex, _, _ = read_h5_scalar_field(np, step, "Ex")
        ey, _, _ = read_h5_scalar_field(np, step, "Ey")
        ez, _, _ = read_h5_scalar_field(np, step, "Ez")

        path_length_s = float(decode_value(attrs.get("path_length_s")))
        interaction_point_s = float(decode_value(attrs.get("interaction_point_s")))
        particle_mean_r = decode_value(attrs.get("particle_mean_r"))
        charge = float(decode_value(attrs.get("particle_total_charge")))

        if not (isinstance(particle_mean_r, list) and len(particle_mean_r) >= 3):
            raise SystemExit("particle_mean_r missing or invalid in HDF5 step.")

        primary_center = (
            float(particle_mean_r[0]),
            float(particle_mean_r[1]),
            float(particle_mean_r[2]),
        )
        centers = [primary_center]
        if snapshot_kind == "active_beambeam":
            interaction_point_beam_z = interaction_point_s - path_length_s
            centers.append(
                (
                    primary_center[0],
                    primary_center[1],
                    2.0 * interaction_point_beam_z - primary_center[2],
                )
            )

        return {
            "step_name": step_name,
            "snapshot_kind": snapshot_kind,
            "origin": origin,
            "spacing": spacing,
            "shape": (rho.shape[2], rho.shape[1], rho.shape[0]),
            "charge": charge,
            "path_length_s": path_length_s,
            "interaction_point_s": interaction_point_s,
            "centers": centers,
            "opalx": {
                "rho": rho,
                "phi": phi,
                "Ex": ex,
                "Ey": ey,
                "Ez": ez,
            },
        }


def compute_error_metrics(np, reference, candidate):
    diff = candidate - reference
    ref_l2 = float(np.sqrt(np.sum(reference * reference)))
    diff_l2 = float(np.sqrt(np.sum(diff * diff)))
    rel_l2 = diff_l2 / ref_l2 if ref_l2 > 0.0 else float("nan")
    max_abs = float(np.max(np.abs(diff)))
    return {
        "max_abs": max_abs,
        "l2": diff_l2,
        "rel_l2": rel_l2,
        "diff": diff,
    }


def compute_volume_integral(np, field, spacing):
    cell_volume = float(spacing[0] * spacing[1] * spacing[2])
    return float(np.sum(field) * cell_volume)


def compute_line_antisymmetry_metrics(np, z_coords, values, symmetry_z):
    mirrored_coords = 2.0 * symmetry_z - z_coords
    valid = (mirrored_coords >= float(np.min(z_coords))) & (mirrored_coords <= float(np.max(z_coords)))
    if not np.any(valid):
        return {
            "max_abs": float("nan"),
            "l2": float("nan"),
            "rel_l2": float("nan"),
            "residual": np.array([], dtype=float),
        }

    mirrored_values = np.interp(mirrored_coords[valid], z_coords, values)
    residual = values[valid] + mirrored_values
    ref_l2 = float(np.sqrt(np.sum(values[valid] * values[valid])))
    diff_l2 = float(np.sqrt(np.sum(residual * residual)))
    return {
        "max_abs": float(np.max(np.abs(residual))),
        "l2": diff_l2,
        "rel_l2": diff_l2 / ref_l2 if ref_l2 > 0.0 else float("nan"),
        "residual": residual,
    }


def central_xz_slice(field, origin, spacing):
    ny = field.shape[1]
    y_index = ny // 2
    projection = field[:, y_index, :]
    x0 = origin[0]
    x1 = origin[0] + field.shape[2] * spacing[0]
    z0 = origin[2]
    z1 = origin[2] + field.shape[0] * spacing[2]
    return projection, (x0, x1, z0, z1)


def central_z_axis_line(np, field, origin, spacing):
    nx = field.shape[2]
    ny = field.shape[1]
    x_index = nx // 2
    y_index = ny // 2
    z_coords = origin[2] + (np.arange(field.shape[0]) + 0.5) * spacing[2]
    values = field[:, y_index, x_index]
    return z_coords, values


def transverse_x_axis_line(np, field, origin, spacing, z_target):
    nx = field.shape[2]
    ny = field.shape[1]
    nz = field.shape[0]
    x_coords = origin[0] + (np.arange(nx) + 0.5) * spacing[0]
    z_coords = origin[2] + (np.arange(nz) + 0.5) * spacing[2]
    y_index = ny // 2
    z_index = int(np.argmin(np.abs(z_coords - z_target)))
    values = field[z_index, y_index, :]
    return x_coords, values


def transverse_y_axis_line(np, field, origin, spacing, z_target):
    nx = field.shape[2]
    ny = field.shape[1]
    nz = field.shape[0]
    y_coords = origin[1] + (np.arange(ny) + 0.5) * spacing[1]
    z_coords = origin[2] + (np.arange(nz) + 0.5) * spacing[2]
    x_index = nx // 2
    z_index = int(np.argmin(np.abs(z_coords - z_target)))
    values = field[z_index, :, x_index]
    return y_coords, values


def derived_output_path(output, suffix):
    if output is None:
        return None
    return output.with_name(f"{output.stem}{suffix}{output.suffix}")


def finalize_figure(plt, fig, output):
    if output is not None:
        fig.savefig(output, dpi=150)
        plt.close(fig)
        return
    plt.show()


def plot_analytic_solution(plt, analytic, origin, spacing, title, output):
    fields = ("rho", "phi")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), squeeze=False, constrained_layout=True)
    fig.suptitle(title)

    for ax, name in zip(axes.flat, fields):
        panel, extent = central_xz_slice(analytic[name], origin, spacing)
        im = ax.imshow(
            panel,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="coolwarm",
        )
        ax.set_title(name)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("z [m]")
        fig.colorbar(im, ax=ax, shrink=0.85)

    finalize_figure(plt, fig, output)


def plot_comparison(
    np,
    plt,
    analytic,
    opalx,
    origin,
    spacing,
    title,
    output,
    components=("rho", "phi"),
    extent_override=None,
    scale_overrides=None,
):
    field_labels = {
        "rho": r"$\rho$",
        "phi": r"$\phi$",
        "Ez": r"$E_z$",
    }
    # components = ("rho", "phi", "Ex", "Ez")
    fig, axes = plt.subplots(
        len(components),
        3,
        figsize=(14, 8),
        squeeze=False,
        constrained_layout=True,
    )
    fig.suptitle(title)

    for row_index, name in enumerate(components):
        metrics = compute_error_metrics(np, analytic[name], opalx[name])
        analytic_panel, extent = central_xz_slice(analytic[name], origin, spacing)
        opalx_panel, _ = central_xz_slice(opalx[name], origin, spacing)
        diff_panel, _ = central_xz_slice(metrics["diff"], origin, spacing)
        if extent_override is not None:
            extent = extent_override

        vmax = max(
            float(abs(analytic_panel).max()),
            float(abs(opalx_panel).max()),
        )
        dmax = float(abs(diff_panel).max())
        if scale_overrides is not None and name in scale_overrides:
            vmax = float(scale_overrides[name]["vmax"])
            dmax = float(scale_overrides[name]["dmax"])

        analytic_ax = axes[row_index][0]
        opalx_ax = axes[row_index][1]
        diff_ax = axes[row_index][2]

        analytic_im = analytic_ax.imshow(
            analytic_panel,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="coolwarm",
            vmin=-vmax if vmax > 0.0 else None,
            vmax=vmax if vmax > 0.0 else None,
        )
        analytic_ax.set_title(f"{field_labels[name]} analytic")
        analytic_ax.set_xlabel("x [m]")
        analytic_ax.set_ylabel("z [m]")

        opalx_im = opalx_ax.imshow(
            opalx_panel,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="coolwarm",
            vmin=-vmax if vmax > 0.0 else None,
            vmax=vmax if vmax > 0.0 else None,
        )
        opalx_ax.set_title(f"{field_labels[name]} OPALX")
        opalx_ax.set_xlabel("x [m]")
        opalx_ax.set_ylabel("z [m]")

        diff_im = diff_ax.imshow(
            diff_panel,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="coolwarm",
            vmin=-dmax if dmax > 0.0 else None,
            vmax=dmax if dmax > 0.0 else None,
        )
        diff_ax.set_title(rf"$\Delta${field_labels[name]}")
        diff_ax.set_xlabel("x [m]")
        diff_ax.set_ylabel("z [m]")

        fig.colorbar(analytic_im, ax=[analytic_ax, opalx_ax], shrink=0.85)
        fig.colorbar(diff_im, ax=diff_ax, shrink=0.85)

        diff_ax.text(
            0.02,
            0.98,
            f"max|Δ|={metrics['max_abs']:.3e}\n"
            f"L2={metrics['l2']:.3e}\n"
            f"rel L2={metrics['rel_l2']:.3e}",
            transform=diff_ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )

    finalize_figure(plt, fig, output)


def plot_rho_z_axis_diagnostic(np, plt, analytic, opalx, origin, spacing, title, output, axis_limits=None):
    z_coords, analytic_line = central_z_axis_line(np, analytic["rho"], origin, spacing)
    _, opalx_line = central_z_axis_line(np, opalx["rho"], origin, spacing)
    delta_line = opalx_line - analytic_line

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), squeeze=False, constrained_layout=True)
    fig.suptitle(title)

    main_ax = axes[0][0]
    main_ax.plot(z_coords, analytic_line, label="analytic", linewidth=2.0)
    main_ax.plot(z_coords, opalx_line, label="OPALX", linewidth=1.8, linestyle="--")
    main_ax.set_xlabel("z [m]")
    main_ax.set_ylabel(r"$\rho$ [C/m$^3$]")
    main_ax.legend()
    main_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        main_ax.set_xlim(*axis_limits["x"])
        main_ax.set_ylim(*axis_limits["main_y"])

    diff_ax = axes[1][0]
    diff_ax.plot(z_coords, delta_line, color="black", linewidth=1.8)
    diff_ax.set_xlabel("z [m]")
    diff_ax.set_ylabel(r"$\Delta\rho$")
    diff_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        diff_ax.set_xlim(*axis_limits["x"])
        diff_ax.set_ylim(*axis_limits["diff_y"])

    finalize_figure(plt, fig, output)


def plot_phi_z_axis_diagnostic(np, plt, analytic, opalx, origin, spacing, title, output, axis_limits=None):
    z_coords, analytic_line = central_z_axis_line(np, analytic["phi"], origin, spacing)
    _, opalx_line = central_z_axis_line(np, opalx["phi"], origin, spacing)
    delta_line = opalx_line - analytic_line

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), squeeze=False, constrained_layout=True)
    fig.suptitle(title)

    main_ax = axes[0][0]
    main_ax.plot(z_coords, analytic_line, label="analytic", linewidth=2.0)
    main_ax.plot(z_coords, opalx_line, label="OPALX", linewidth=1.8, linestyle="--")
    main_ax.set_xlabel("z [m]")
    main_ax.set_ylabel(r"$\phi$ [V]")
    main_ax.legend()
    main_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        main_ax.set_xlim(*axis_limits["x"])
        main_ax.set_ylim(*axis_limits["main_y"])

    diff_ax = axes[1][0]
    diff_ax.plot(z_coords, delta_line, color="black", linewidth=1.8)
    diff_ax.set_xlabel("z [m]")
    diff_ax.set_ylabel(r"$\Delta\phi$")
    diff_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        diff_ax.set_xlim(*axis_limits["x"])
        diff_ax.set_ylim(*axis_limits["diff_y"])

    finalize_figure(plt, fig, output)


def plot_ez_z_axis_diagnostic(np, plt, analytic, opalx, origin, spacing, title, output, axis_limits=None):
    z_coords, analytic_line = central_z_axis_line(np, analytic["Ez"], origin, spacing)
    _, opalx_line = central_z_axis_line(np, opalx["Ez"], origin, spacing)
    delta_line = opalx_line - analytic_line

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), squeeze=False, constrained_layout=True)
    fig.suptitle(title)

    main_ax = axes[0][0]
    main_ax.plot(z_coords, analytic_line, label="analytic", linewidth=2.0)
    main_ax.plot(z_coords, opalx_line, label="OPALX", linewidth=1.8, linestyle="--")
    main_ax.set_xlabel("z [m]")
    main_ax.set_ylabel(r"$E_z$ [V/m]")
    main_ax.legend()
    main_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        main_ax.set_xlim(*axis_limits["x"])
        main_ax.set_ylim(*axis_limits["main_y"])

    diff_ax = axes[1][0]
    diff_ax.plot(z_coords, delta_line, color="black", linewidth=1.8)
    diff_ax.set_xlabel("z [m]")
    diff_ax.set_ylabel(r"$\Delta E_z$")
    diff_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        diff_ax.set_xlim(*axis_limits["x"])
        diff_ax.set_ylim(*axis_limits["diff_y"])

    finalize_figure(plt, fig, output)


def plot_ex_x_axis_diagnostic(
    np, plt, analytic, opalx, origin, spacing, title, output, z_target, axis_limits=None
):
    x_coords, analytic_line = transverse_x_axis_line(np, analytic["Ex"], origin, spacing, z_target)
    _, opalx_line = transverse_x_axis_line(np, opalx["Ex"], origin, spacing, z_target)
    delta_line = opalx_line - analytic_line

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), squeeze=False, constrained_layout=True)
    fig.suptitle(title)

    main_ax = axes[0][0]
    main_ax.plot(x_coords, analytic_line, label="analytic", linewidth=2.0)
    main_ax.plot(x_coords, opalx_line, label="OPALX", linewidth=1.8, linestyle="--")
    main_ax.set_xlabel("x [m]")
    main_ax.set_ylabel(r"$E_x$ [V/m]")
    main_ax.legend()
    main_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        main_ax.set_xlim(*axis_limits["x"])
        main_ax.set_ylim(*axis_limits["main_y"])

    diff_ax = axes[1][0]
    diff_ax.plot(x_coords, delta_line, color="black", linewidth=1.8)
    diff_ax.set_xlabel("x [m]")
    diff_ax.set_ylabel(r"$\Delta E_x$")
    diff_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        diff_ax.set_xlim(*axis_limits["x"])
        diff_ax.set_ylim(*axis_limits["diff_y"])

    finalize_figure(plt, fig, output)


def plot_ey_y_axis_diagnostic(
    np, plt, analytic, opalx, origin, spacing, title, output, z_target, axis_limits=None
):
    y_coords, analytic_line = transverse_y_axis_line(np, analytic["Ey"], origin, spacing, z_target)
    _, opalx_line = transverse_y_axis_line(np, opalx["Ey"], origin, spacing, z_target)
    delta_line = opalx_line - analytic_line

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), squeeze=False, constrained_layout=True)
    fig.suptitle(title)

    main_ax = axes[0][0]
    main_ax.plot(y_coords, analytic_line, label="analytic", linewidth=2.0)
    main_ax.plot(y_coords, opalx_line, label="OPALX", linewidth=1.8, linestyle="--")
    main_ax.set_xlabel("y [m]")
    main_ax.set_ylabel(r"$E_y$ [V/m]")
    main_ax.legend()
    main_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        main_ax.set_xlim(*axis_limits["x"])
        main_ax.set_ylim(*axis_limits["main_y"])

    diff_ax = axes[1][0]
    diff_ax.plot(y_coords, delta_line, color="black", linewidth=1.8)
    diff_ax.set_xlabel("y [m]")
    diff_ax.set_ylabel(r"$\Delta E_y$")
    diff_ax.grid(True, alpha=0.3)
    if axis_limits is not None:
        diff_ax.set_xlim(*axis_limits["x"])
        diff_ax.set_ylim(*axis_limits["diff_y"])

    finalize_figure(plt, fig, output)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Manufactured-solution tool for symmetric 3D Gaussian BeamBeam test cases."
    )
    parser.add_argument("--sigma", type=float, default=1.0e-3, help="Gaussian sigma [m]")
    parser.add_argument(
        "--charge",
        type=float,
        default=-1.0e-15,
        help="Charge per bunch [C] for pure analytic mode",
    )
    parser.add_argument(
        "--half-separation",
        type=float,
        default=2.0e-3,
        help="Half distance from IP to each bunch center [m] in analytic mode",
    )
    parser.add_argument("--nx", type=int, default=64, help="Grid points in x for analytic mode")
    parser.add_argument("--ny", type=int, default=64, help="Grid points in y for analytic mode")
    parser.add_argument("--nz", type=int, default=64, help="Grid points in z for analytic mode")
    parser.add_argument("--x-span", type=float, help="Total x span [m] for analytic mode")
    parser.add_argument("--y-span", type=float, help="Total y span [m] for analytic mode")
    parser.add_argument("--z-span", type=float, help="Total z span [m] for analytic mode")
    parser.add_argument(
        "--compare-h5",
        type=Path,
        help="Compare against one OPALX HDF5 BeamBeam diagnostics file (OPEN solver case)",
    )
    parser.add_argument("--compare-rho-dump", type=Path, help="Compare an OPALX ASCII rho dump.")
    parser.add_argument("--compare-phi-dump", type=Path, help="Companion OPALX ASCII phi dump.")
    parser.add_argument("--compare-e-dump", type=Path, help="Companion OPALX ASCII E-vector dump.")
    parser.add_argument("--step", type=int, help="Global step to read from --compare-h5")
    parser.add_argument(
        "--state",
        help="Optional raw HDF5 state selector, e.g. active_beambeam or normal_tracking",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional PNG output path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    np = load_numpy()

    if args.compare_rho_dump is not None:
        if args.compare_phi_dump is None or args.compare_e_dump is None:
            raise SystemExit("--compare-rho-dump requires --compare-phi-dump and --compare-e-dump")

        setup = manufactured_setup_from_ascii(
            np, args.compare_rho_dump, args.compare_phi_dump, args.compare_e_dump)
        origin = setup["origin"]
        spacing = setup["spacing"]
        shape = setup["shape"]
        charge = setup["charge"]
        centers = setup["centers"]

        X, Y, Z = build_grid(np, origin, spacing, shape)
        analytic = gaussian_pair_fields(np, X, Y, Z, args.sigma, charge, centers)

        print(f"rho dump: {args.compare_rho_dump}")
        print(f"phi dump: {args.compare_phi_dump}")
        print(f"E dump: {args.compare_e_dump}")
        print(f"state: {setup['snapshot_kind']}")
        print(f"charge per bunch: {charge:.6e} C")
        for index, center in enumerate(centers, start=1):
            print(
                f"center{index} (beam-local frame): "
                f"({center[0]:.6e}, {center[1]:.6e}, {center[2]:.6e})"
            )

        for name in ("rho", "phi", "Ex", "Ey", "Ez"):
            metrics = compute_error_metrics(np, analytic[name], setup["opalx"][name])
            print(
                f"{name}: max|Δ|={metrics['max_abs']:.6e}, "
                f"L2={metrics['l2']:.6e}, relL2={metrics['rel_l2']:.6e}"
            )

        ip_local_z = setup["interaction_point_local_z"]
        z_coords, analytic_ez_line = central_z_axis_line(np, analytic["Ez"], origin, spacing)
        _, opalx_ez_line = central_z_axis_line(np, setup["opalx"]["Ez"], origin, spacing)
        ip_index = int(np.argmin(np.abs(z_coords - ip_local_z)))
        print(
            f"Ez(nearest grid sample to IP) on axis: analytic={analytic_ez_line[ip_index]:.6e}, "
            f"OPALX={opalx_ez_line[ip_index]:.6e} [V/m]"
        )
        return 0

    if args.compare_h5 is not None:
        plt = load_matplotlib()
        if args.step is None:
            raise SystemExit("--step is required with --compare-h5")

        setup = manufactured_setup_from_h5(np, args.compare_h5, args.step, args.state)
        origin = setup["origin"]
        spacing = setup["spacing"]
        shape = setup["shape"]
        charge = setup["charge"]
        centers = setup["centers"]

        X, Y, Z = build_grid(np, origin, spacing, shape)
        analytic = gaussian_pair_fields(np, X, Y, Z, args.sigma, charge, centers)

        print(f"file: {args.compare_h5}")
        print(f"step: {args.step}")
        print(f"state: {setup['snapshot_kind']}")
        print(f"charge per bunch: {charge:.6e} C")
        for index, center in enumerate(centers, start=1):
            print(
                f"center{index} (beam frame): "
                f"({center[0]:.6e}, {center[1]:.6e}, {center[2]:.6e})"
            )

        for name in ("rho", "phi", "Ex", "Ey", "Ez"):
            metrics = compute_error_metrics(np, analytic[name], setup["opalx"][name])
            print(
                f"{name}: max|Δ|={metrics['max_abs']:.6e}, "
                f"L2={metrics['l2']:.6e}, relL2={metrics['rel_l2']:.6e}"
            )

        ip_beam_z = setup["interaction_point_s"] - setup["path_length_s"]
        z_coords, analytic_ez_line = central_z_axis_line(np, analytic["Ez"], origin, spacing)
        _, opalx_ez_line = central_z_axis_line(np, setup["opalx"]["Ez"], origin, spacing)
        analytic_integral = compute_volume_integral(np, analytic["Ez"], spacing)
        opalx_integral = compute_volume_integral(np, setup["opalx"]["Ez"], spacing)
        analytic_antisym = compute_line_antisymmetry_metrics(np, z_coords, analytic_ez_line, ip_beam_z)
        opalx_antisym = compute_line_antisymmetry_metrics(np, z_coords, opalx_ez_line, ip_beam_z)
        ip_index = int(np.argmin(np.abs(z_coords - ip_beam_z)))
        print(
            f"Ez integral: analytic={analytic_integral:.6e}, "
            f"OPALX={opalx_integral:.6e} [V m^2]"
        )
        print(
            f"Ez(IP) on axis: analytic={analytic_ez_line[ip_index]:.6e}, "
            f"OPALX={opalx_ez_line[ip_index]:.6e} [V/m]"
        )
        print(
            f"Ez antisymmetry on axis: analytic max|res|={analytic_antisym['max_abs']:.6e}, "
            f"OPALX max|res|={opalx_antisym['max_abs']:.6e}"
        )

        x_coords, analytic_ex_line = transverse_x_axis_line(np, analytic["Ex"], origin, spacing, ip_beam_z)
        _, opalx_ex_line = transverse_x_axis_line(np, setup["opalx"]["Ex"], origin, spacing, ip_beam_z)
        analytic_ex_antisym = compute_line_antisymmetry_metrics(np, x_coords, analytic_ex_line, 0.0)
        opalx_ex_antisym = compute_line_antisymmetry_metrics(np, x_coords, opalx_ex_line, 0.0)
        x0_index = int(np.argmin(np.abs(x_coords)))
        print(
            f"Ex(x, z≈IP) on axis: analytic={analytic_ex_line[x0_index]:.6e}, "
            f"OPALX={opalx_ex_line[x0_index]:.6e} [V/m]"
        )
        print(
            f"Ex antisymmetry at IP plane: analytic max|res|={analytic_ex_antisym['max_abs']:.6e}, "
            f"OPALX max|res|={opalx_ex_antisym['max_abs']:.6e}"
        )

        y_coords, analytic_ey_line = transverse_y_axis_line(np, analytic["Ey"], origin, spacing, ip_beam_z)
        _, opalx_ey_line = transverse_y_axis_line(np, setup["opalx"]["Ey"], origin, spacing, ip_beam_z)
        analytic_ey_antisym = compute_line_antisymmetry_metrics(np, y_coords, analytic_ey_line, 0.0)
        opalx_ey_antisym = compute_line_antisymmetry_metrics(np, y_coords, opalx_ey_line, 0.0)
        y0_index = int(np.argmin(np.abs(y_coords)))
        print(
            f"Ey(y, z≈IP) on axis: analytic={analytic_ey_line[y0_index]:.6e}, "
            f"OPALX={opalx_ey_line[y0_index]:.6e} [V/m]"
        )
        print(
            f"Ey antisymmetry at IP plane: analytic max|res|={analytic_ey_antisym['max_abs']:.6e}, "
            f"OPALX max|res|={opalx_ey_antisym['max_abs']:.6e}"
        )

        title = (
            f"Manufactured solution vs OPALX | "
            f"step {args.step} | {setup['snapshot_kind']} | "
            f"$Q$={charge:.2e} C, $\\sigma$={args.sigma:.2e} m"
        )
        plot_comparison(np, plt, analytic, setup["opalx"], origin, spacing, title, args.output)
        rho_line_output = derived_output_path(args.output, "-rho-z-axis")
        plot_rho_z_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            f"{title}",
            rho_line_output,
        )
        phi_line_output = derived_output_path(args.output, "-phi-z-axis")
        plot_phi_z_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            f"{title}",
            phi_line_output,
        )
        ez_line_output = derived_output_path(args.output, "-ez-z-axis")
        plot_ez_z_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            f"{title}",
            ez_line_output,
        )
        ex_line_output = derived_output_path(args.output, "-ex-x-axis")
        plot_ex_x_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            f"{title}",
            ex_line_output,
            ip_beam_z,
        )
        ey_line_output = derived_output_path(args.output, "-ey-y-axis")
        plot_ey_y_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            f"{title}",
            ey_line_output,
            ip_beam_z,
        )
        return 0

    plt = load_matplotlib()
    origin, spacing, shape = default_grid(args)
    X, Y, Z = build_grid(np, origin, spacing, shape)
    centers = [
        (0.0, 0.0, -args.half_separation),
        (0.0, 0.0, +args.half_separation),
    ]
    analytic = gaussian_pair_fields(np, X, Y, Z, args.sigma, args.charge, centers)

    print("analytic symmetric Gaussian pair")
    print(f"sigma = {args.sigma:.6e} m")
    print(f"charge per bunch = {args.charge:.6e} C")
    print(f"centers = {centers}")
    print(f"origin = {origin}")
    print(f"spacing = {spacing}")
    print(f"shape = {shape}")

    title = (
        "BeamBeam manufactured solution: symmetric 3D Gaussian pair"
        f" | $Q$={args.charge:.3e} C, $\\sigma$={args.sigma:.3e} m"
    )
    plot_analytic_solution(plt, analytic, origin, spacing, title, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
