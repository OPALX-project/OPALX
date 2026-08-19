#!/usr/bin/env python3
"""Evaluate E and B for two rigid, counter-propagating Gaussian sources.

The interaction point is the origin. At each requested snapshot, two identical
primary electron-beam centroids are placed at z=-d/2 and z=+d/2 and move toward
the IP with velocities +beta*c and -beta*c, respectively. Here d is expressed
in units of the lab-frame longitudinal sigma_z.

Each source is evaluated in its rest frame as a triaxial Gaussian and Lorentz
transformed back to the lab. Pair particles are passive probes and do not
contribute to these fields.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
ROOT = MODEL_DIR.parents[1]
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_CONFIG = MODEL_DIR / "parameters.json"
DEFAULT_OUTPUT_DIR = MODEL_DIR / "outputs" / "fields"


@dataclass(frozen=True)
class ModelParameters:
    kinetic_energy_MeV: float
    electrons_per_bunch: float
    sigma_x_m: float
    sigma_y_m: float
    sigma_z_m: float
    centroid_separations_sigma_z: tuple[float, ...]
    plane_y_m: float
    x_half_width_sigma_x: float
    z_half_width_sigma_z: float
    nx: int
    nz: int
    quadrature_chunk_size: int


def load_track12_module():
    """Load the established scalar boosted-Gaussian evaluator for cross-checks."""
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_parameters(path: Path) -> ModelParameters:
    with path.open("r", encoding="utf-8") as stream:
        config = json.load(stream)

    source = config["source"]
    snapshots = config["snapshots"]
    grid = config["grid"]
    parameters = ModelParameters(
        kinetic_energy_MeV=float(source["kinetic_energy_MeV"]),
        electrons_per_bunch=float(source["electrons_per_bunch"]),
        sigma_x_m=float(source["sigma_x_m"]),
        sigma_y_m=float(source["sigma_y_m"]),
        sigma_z_m=float(source["sigma_z_m"]),
        centroid_separations_sigma_z=tuple(
            float(value) for value in snapshots["centroid_separations_sigma_z"]
        ),
        plane_y_m=float(snapshots["plane_y_m"]),
        x_half_width_sigma_x=float(grid["x_half_width_sigma_x"]),
        z_half_width_sigma_z=float(grid["z_half_width_sigma_z"]),
        nx=int(grid["nx"]),
        nz=int(grid["nz"]),
        quadrature_chunk_size=int(grid["quadrature_chunk_size"]),
    )
    validate_parameters(parameters)
    return parameters


def validate_parameters(parameters: ModelParameters) -> None:
    positive_values = {
        "kinetic_energy_MeV": parameters.kinetic_energy_MeV,
        "electrons_per_bunch": parameters.electrons_per_bunch,
        "sigma_x_m": parameters.sigma_x_m,
        "sigma_y_m": parameters.sigma_y_m,
        "sigma_z_m": parameters.sigma_z_m,
        "x_half_width_sigma_x": parameters.x_half_width_sigma_x,
        "z_half_width_sigma_z": parameters.z_half_width_sigma_z,
        "quadrature_chunk_size": parameters.quadrature_chunk_size,
    }
    for name, value in positive_values.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive, got {value}")
    if parameters.nx < 3 or parameters.nz < 3:
        raise ValueError("nx and nz must both be at least 3")
    if parameters.nx % 2 == 0 or parameters.nz % 2 == 0:
        raise ValueError("nx and nz must be odd so that the IP is sampled exactly")
    if not parameters.centroid_separations_sigma_z:
        raise ValueError("at least one centroid separation is required")
    if any(value < 0 for value in parameters.centroid_separations_sigma_z):
        raise ValueError("centroid separations must be non-negative")


def make_sources(track12, parameters: ModelParameters, separation_sigma_z: float):
    beta = track12.beta_from_kinetic_energy(parameters.kinetic_energy_MeV)
    charge_C = -parameters.electrons_per_bunch * track12.ELEMENTARY_CHARGE_C
    half_separation_m = 0.5 * separation_sigma_z * parameters.sigma_z_m
    sigma_lab_m = (parameters.sigma_x_m, parameters.sigma_y_m, parameters.sigma_z_m)

    left_to_right = track12.RigidAnisotropicGaussianSource(
        "left_to_right",
        charge_C=charge_C,
        sigma_lab_m=sigma_lab_m,
        beta_z=+beta,
        center_t0_m=np.array([0.0, 0.0, -half_separation_m]),
    )
    right_to_left = track12.RigidAnisotropicGaussianSource(
        "right_to_left",
        charge_C=charge_C,
        sigma_lab_m=sigma_lab_m,
        beta_z=-beta,
        center_t0_m=np.array([0.0, 0.0, +half_separation_m]),
    )
    return left_to_right, right_to_left


def anisotropic_gaussian_rest_fields_batch(
    positions_rest_m: np.ndarray,
    charge_C: float,
    sigma_rest_m: np.ndarray,
    *,
    epsilon_0: float,
    log_quad_nodes: np.ndarray,
    chunk_size: int,
) -> np.ndarray:
    """Vectorized form of the established one-dimensional Gaussian quadrature."""
    positions_rest_m = np.asarray(positions_rest_m, dtype=float)
    sigma2 = np.asarray(sigma_rest_m, dtype=float) ** 2
    scale = float(np.exp(np.mean(np.log(sigma2))))
    t = scale * np.exp(log_quad_nodes)
    denominators = sigma2[None, :] + t[:, None]
    denominator_product = np.sqrt(np.prod(denominators, axis=1))
    integration_kernel = (
        t[:, None] / (denominator_product[:, None] * denominators)
    )
    coefficient = charge_C / (4.0 * np.pi * epsilon_0 * np.sqrt(2.0 * np.pi))
    fields = np.empty_like(positions_rest_m)

    for begin in range(0, len(positions_rest_m), chunk_size):
        end = min(begin + chunk_size, len(positions_rest_m))
        positions = positions_rest_m[begin:end]
        exponent = -0.5 * np.sum(
            positions[:, None, :] ** 2 / denominators[None, :, :], axis=2
        )
        integrals = np.trapezoid(
            np.exp(exponent)[:, :, None] * integration_kernel[None, :, :],
            log_quad_nodes,
            axis=1,
        )
        fields[begin:end] = coefficient * positions * integrals
    return fields


def source_lab_fields_batch(
    positions_m: np.ndarray,
    source,
    track12,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    displacement = np.asarray(positions_m, dtype=float) - source.center(0.0)[None, :]
    positions_rest = displacement.copy()
    positions_rest[:, 2] *= source.gamma
    e_rest = anisotropic_gaussian_rest_fields_batch(
        positions_rest,
        source.charge_C,
        source.sigma_rest,
        epsilon_0=track12.EPSILON_0,
        log_quad_nodes=track12.LOG_QUAD_NODES,
        chunk_size=chunk_size,
    )
    e_lab = e_rest.copy()
    e_lab[:, :2] *= source.gamma
    velocity = np.broadcast_to(source.beta * track12.C_LIGHT, e_lab.shape)
    b_lab = np.cross(velocity, e_lab) / track12.C_LIGHT**2
    return e_lab, b_lab


def total_lab_fields_batch(
    positions_m: np.ndarray,
    sources,
    track12,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    e_total = np.zeros_like(positions_m, dtype=float)
    b_total = np.zeros_like(positions_m, dtype=float)
    for source in sources:
        e_field, b_field = source_lab_fields_batch(
            positions_m, source, track12, chunk_size
        )
        e_total += e_field
        b_total += b_field
    return e_total, b_total


def scalar_cross_check(
    positions_m: np.ndarray,
    e_batch: np.ndarray,
    b_batch: np.ndarray,
    sources,
    track12,
) -> dict[str, float]:
    sample_indices = np.unique(
        np.linspace(0, len(positions_m) - 1, min(17, len(positions_m)), dtype=int)
    )
    e_scalar = []
    b_scalar = []
    for index in sample_indices:
        e_field, b_field = track12.anisotropic_total_lab_fields(
            positions_m[index], 0.0, sources
        )
        e_scalar.append(e_field)
        b_scalar.append(b_field)
    e_scalar = np.asarray(e_scalar)
    b_scalar = np.asarray(b_scalar)
    e_error = e_batch[sample_indices] - e_scalar
    b_error = b_batch[sample_indices] - b_scalar
    e_scale = max(float(np.max(np.abs(e_scalar))), np.finfo(float).tiny)
    b_scale = max(float(np.max(np.abs(b_scalar))), np.finfo(float).tiny)
    return {
        "sample_count": int(len(sample_indices)),
        "max_abs_E_error_V_per_m": float(np.max(np.abs(e_error))),
        "max_rel_E_error": float(np.max(np.abs(e_error)) / e_scale),
        "max_abs_B_error_T": float(np.max(np.abs(b_error))),
        "max_rel_B_error": float(np.max(np.abs(b_error)) / b_scale),
    }


def separation_label(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p") + "sigma"


def calculate_snapshot(track12, parameters: ModelParameters, separation_sigma_z: float):
    x_values = np.linspace(
        -parameters.x_half_width_sigma_x * parameters.sigma_x_m,
        +parameters.x_half_width_sigma_x * parameters.sigma_x_m,
        parameters.nx,
    )
    z_values = np.linspace(
        -parameters.z_half_width_sigma_z * parameters.sigma_z_m,
        +parameters.z_half_width_sigma_z * parameters.sigma_z_m,
        parameters.nz,
    )
    z_grid, x_grid = np.meshgrid(z_values, x_values, indexing="ij")
    positions = np.column_stack(
        (
            x_grid.ravel(),
            np.full(x_grid.size, parameters.plane_y_m),
            z_grid.ravel(),
        )
    )
    sources = make_sources(track12, parameters, separation_sigma_z)
    e_field, b_field = total_lab_fields_batch(
        positions, sources, track12, parameters.quadrature_chunk_size
    )
    data = pd.DataFrame(
        {
            "separation_sigma_z": separation_sigma_z,
            "separation_m": separation_sigma_z * parameters.sigma_z_m,
            "x_m": positions[:, 0],
            "y_m": positions[:, 1],
            "z_minus_ip_m": positions[:, 2],
            "Ex_V_per_m": e_field[:, 0],
            "Ey_V_per_m": e_field[:, 1],
            "Ez_V_per_m": e_field[:, 2],
            "Bx_T": b_field[:, 0],
            "By_T": b_field[:, 1],
            "Bz_T": b_field[:, 2],
            "E_abs_V_per_m": np.linalg.norm(e_field, axis=1),
            "B_abs_T": np.linalg.norm(b_field, axis=1),
        }
    )
    shape = (parameters.nz, parameters.nx)
    cross_check = scalar_cross_check(positions, e_field, b_field, sources, track12)
    return data, e_field.reshape(*shape, 3), b_field.reshape(*shape, 3), cross_check


def snapshot_summary(
    separation_sigma_z: float,
    data: pd.DataFrame,
    e_grid: np.ndarray,
    b_grid: np.ndarray,
    cross_check: dict[str, float],
) -> dict[str, float]:
    e_scale = max(float(np.max(np.abs(e_grid))), np.finfo(float).tiny)
    b_scale = max(float(np.max(np.abs(b_grid))), np.finfo(float).tiny)
    center_z = e_grid.shape[0] // 2
    center_x = e_grid.shape[1] // 2
    return {
        "separation_sigma_z": separation_sigma_z,
        "separation_m": float(data["separation_m"].iloc[0]),
        "max_E_abs_V_per_m": float(data["E_abs_V_per_m"].max()),
        "max_B_abs_T": float(data["B_abs_T"].max()),
        "ip_E_abs_V_per_m": float(np.linalg.norm(e_grid[center_z, center_x])),
        "ip_B_abs_T": float(np.linalg.norm(b_grid[center_z, center_x])),
        "E_inversion_odd_residual": float(
            np.max(np.abs(e_grid + e_grid[::-1, ::-1, :])) / e_scale
        ),
        "B_inversion_even_residual": float(
            np.max(np.abs(b_grid - b_grid[::-1, ::-1, :])) / b_scale
        ),
        **cross_check,
    }


def engineering_scale(maximum: float) -> tuple[float, str]:
    if maximum >= 1.0e12:
        return 1.0e12, "T"
    if maximum >= 1.0e9:
        return 1.0e9, "G"
    if maximum >= 1.0e6:
        return 1.0e6, "M"
    if maximum >= 1.0e3:
        return 1.0e3, "k"
    if maximum >= 1.0:
        return 1.0, ""
    if maximum >= 1.0e-3:
        return 1.0e-3, "m"
    if maximum >= 1.0e-6:
        return 1.0e-6, "u"
    return 1.0, ""


def configure_matplotlib(output_dir: Path) -> None:
    cache_dir = output_dir / ".plot-cache" / "matplotlib"
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    cache_dir.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def plot_snapshots(snapshots, parameters: ModelParameters, output_dir: Path) -> Path:
    configure_matplotlib(output_dir)
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    components = (
        ("Ex", 0, "E", "V/m"),
        ("Ez", 2, "E", "V/m"),
        ("By", 1, "B", "T"),
    )
    maxima = {}
    for label, axis, field_kind, _unit in components:
        maxima[label] = max(
            float(np.max(np.abs(e_grid[..., axis] if field_kind == "E" else b_grid[..., axis])))
            for _separation, _data, e_grid, b_grid, _check in snapshots
        )
    scales = {label: engineering_scale(maximum) for label, maximum in maxima.items()}

    nrows = len(snapshots)
    fig, axes = plt.subplots(
        nrows,
        len(components),
        figsize=(11.0, 2.55 * nrows + 0.8),
        sharex=True,
        sharey=True,
        constrained_layout=True,
        squeeze=False,
    )
    x_normalized = np.linspace(
        -parameters.x_half_width_sigma_x, parameters.x_half_width_sigma_x, parameters.nx
    )
    z_normalized = np.linspace(
        -parameters.z_half_width_sigma_z, parameters.z_half_width_sigma_z, parameters.nz
    )
    images = []
    for row, (separation, _data, e_grid, b_grid, _check) in enumerate(snapshots):
        for column, (label, axis, field_kind, unit) in enumerate(components):
            field = e_grid[..., axis] if field_kind == "E" else b_grid[..., axis]
            scale, prefix = scales[label]
            limit = maxima[label] / scale
            if limit == 0.0:
                limit = 1.0
            image = axes[row, column].pcolormesh(
                z_normalized,
                x_normalized,
                (field / scale).T,
                shading="auto",
                cmap="RdBu_r",
                norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
                rasterized=True,
            )
            if row == 0:
                axes[row, column].set_title(f"${label}$ [{prefix}{unit}]")
                images.append(image)
            half_separation = 0.5 * separation
            axes[row, column].plot(
                [-half_separation, +half_separation],
                [0.0, 0.0],
                marker="x",
                linestyle="none",
                color="black",
                markersize=5,
                markeredgewidth=1.1,
            )
            axes[row, column].axvline(0.0, color="black", linewidth=0.55, alpha=0.45)
            axes[row, column].grid(False)
        axes[row, 0].set_ylabel(
            f"$x/\\sigma_x$\n$d={separation:g}\\sigma_z$"
        )

    for column, image in enumerate(images):
        fig.colorbar(image, ax=axes[:, column], shrink=0.92, pad=0.02)
        axes[-1, column].set_xlabel("$(z-z_\\mathrm{IP})/\\sigma_z$")
    fig.suptitle(
        "Rigid counter-propagating primary beams: fields in the $y=0$ plane\n"
        "black crosses: primary-beam centroids",
        fontsize=12,
    )
    output_path = output_dir / "rigid_two_gaussian_fields_xz.png"
    fig.savefig(output_path, dpi=240)
    plt.close(fig)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    parameters = load_parameters(args.config)
    track12 = load_track12_module()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    snapshots = []
    summaries = []
    for separation in parameters.centroid_separations_sigma_z:
        data, e_grid, b_grid, cross_check = calculate_snapshot(
            track12, parameters, separation
        )
        data_path = args.output_dir / f"fields_{separation_label(separation)}.csv"
        data.to_csv(data_path, index=False)
        snapshots.append((separation, data, e_grid, b_grid, cross_check))
        summaries.append(snapshot_summary(separation, data, e_grid, b_grid, cross_check))
        print(f"wrote {data_path}")

    summary_frame = pd.DataFrame(summaries)
    summary_csv = args.output_dir / "field_snapshot_summary.csv"
    summary_frame.to_csv(summary_csv, index=False)

    beta = track12.beta_from_kinetic_energy(parameters.kinetic_energy_MeV)
    run_summary = {
        "config": str(args.config.resolve()),
        "parameters": asdict(parameters),
        "derived": {
            "beta": beta,
            "gamma": 1.0 / math.sqrt(1.0 - beta * beta),
            "charge_per_bunch_C": (
                -parameters.electrons_per_bunch * track12.ELEMENTARY_CHARGE_C
            ),
        },
        "model": {
            "interaction_point_m": [0.0, 0.0, 0.0],
            "centers_m": "(0,0,-d/2) moving +z and (0,0,+d/2) moving -z",
            "sigma_z_definition": "lab-frame rms longitudinal size",
            "pair_fields_included": False,
            "retardation": "exact uniform-motion Lorentz transform of each rigid source",
        },
        "products": [
            f"fields_{separation_label(value)}.csv"
            for value in parameters.centroid_separations_sigma_z
        ],
    }
    with (args.output_dir / "run_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(run_summary, stream, indent=2)
        stream.write("\n")

    plot_path = None
    if not args.no_plots:
        plot_path = plot_snapshots(snapshots, parameters, args.output_dir)

    print(summary_frame.to_string(index=False, float_format=lambda value: f"{value:.6e}"))
    print(f"wrote {summary_csv}")
    if plot_path is not None:
        print(f"wrote {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
