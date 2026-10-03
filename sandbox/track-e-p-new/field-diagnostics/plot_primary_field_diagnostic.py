#!/usr/bin/env python3
"""Visualize a primary-only OPALX field dump against the manufactured field."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MODEL_SCRIPT = ROOT / "sandbox/reduced-order-model/python/rigid_two_gaussian_fields.py"
COMPARE_SCRIPT = ROOT / "sandbox/reduced-order-model/python/compare_opalx_fields.py"
DEFAULT_CASE = HERE / "outputs/0sigma"
DEFAULT_CONFIG = HERE / "primary_field_parameters.json"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_probe_positions(path: Path) -> np.ndarray:
    frame = pd.read_csv(path, sep=r"\s+", skiprows=1)
    return frame[["x", "y", "z"]].to_numpy(dtype=float)


def field_magnitude(group: h5py.Group, names: tuple[str, str, str]) -> tuple[np.ndarray, np.ndarray]:
    values = np.column_stack([np.asarray(group[name], dtype=float) for name in names])
    finite = np.isfinite(values).all(axis=1)
    magnitude = np.full(len(values), np.nan)
    magnitude[finite] = np.linalg.norm(values[finite], axis=1)
    return magnitude, finite


def positive_limits(values: np.ndarray) -> tuple[float, float]:
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        return 1.0, 10.0
    vmax = float(np.max(positive))
    vmin = max(float(np.quantile(positive, 0.01)), vmax * 1.0e-7)
    return vmin, vmax


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()

    field_model = load_module("bb_field_model", MODEL_SCRIPT)
    compare = load_module("bb_field_compare", COMPARE_SCRIPT)
    parameters = field_model.load_parameters(args.config)
    track12 = field_model.load_track12_module()

    h5_path = next(args.case_dir.glob("*_c1.h5"))
    stat_path = next(args.case_dir.glob("*_c0.stat"))
    probe_path = next(args.case_dir.glob("*_probes.fromfile"))
    input_path = next(args.case_dir.glob("*.in"))
    deck = input_path.read_text(encoding="utf-8")
    green_match = re.search(r"GREENSF\s*=\s*([A-Za-z]+)", deck)
    nxy_match = re.search(r"REAL\s+NXY\s*=\s*(\d+)", deck)
    nz_match = re.search(r"REAL\s+NZ\s*=\s*(\d+)", deck)
    if green_match is None or nxy_match is None or nz_match is None:
        raise ValueError(f"cannot read Green function and mesh from {input_path}")
    green_function = green_match.group(1).upper()
    mesh = [int(nxy_match.group(1)), int(nxy_match.group(1)), int(nz_match.group(1))]
    positions_local = load_probe_positions(probe_path)
    expected = parameters.nx * parameters.nz
    if len(positions_local) != expected:
        raise ValueError(f"expected {expected} probes, found {len(positions_local)}")

    ip_s_m = 0.008
    source_centroid_abs_m = compare.read_source_centroid_from_stat(stat_path, 0)
    source_center_from_ip_m = source_centroid_abs_m - ip_s_m
    positions_from_ip = positions_local.copy()

    beta = track12.beta_from_kinetic_energy(parameters.kinetic_energy_MeV)
    source = track12.RigidAnisotropicGaussianSource(
        "physical",
        charge_C=-parameters.electrons_per_bunch * track12.ELEMENTARY_CHARGE_C,
        sigma_lab_m=(parameters.sigma_x_m, parameters.sigma_y_m, parameters.sigma_z_m),
        beta_z=+beta,
        center_t0_m=np.array([0.0, 0.0, source_center_from_ip_m]),
        cutoff_sigma=3.0,
    )
    e_manufactured, b_manufactured = field_model.total_lab_fields_batch(
        positions_from_ip, (source,), track12, parameters.quadrature_chunk_size
    )
    e_manufactured_abs = np.linalg.norm(e_manufactured, axis=1)
    b_manufactured_abs = np.linalg.norm(b_manufactured, axis=1)

    with h5py.File(h5_path, "r") as h5_file:
        group = h5_file["Step#0"]
        e_opalx_abs, e_finite = field_magnitude(group, ("Ex", "Ey", "Ez"))
        b_opalx_abs, b_finite = field_magnitude(group, ("Bx", "By", "Bz"))
        opalx_components = {
            name: np.asarray(group[name], dtype=float)
            for name in ("Ex", "Ey", "Ez", "Bx", "By", "Bz")
        }

    e_opalx = np.column_stack(
        [opalx_components[name] for name in ("Ex", "Ey", "Ez")]
    )
    b_opalx = np.column_stack(
        [opalx_components[name] for name in ("Bx", "By", "Bz")]
    )

    table = pd.DataFrame(
        {
            "x_m": positions_from_ip[:, 0],
            "y_m": positions_from_ip[:, 1],
            "z_minus_ip_m": positions_from_ip[:, 2],
            **{f"{name}_opalx": values for name, values in opalx_components.items()},
            **{
                f"{name}_manufactured": values
                for name, values in zip(("Ex", "Ey", "Ez"), e_manufactured.T, strict=True)
            },
            **{
                f"{name}_manufactured": values
                for name, values in zip(("Bx", "By", "Bz"), b_manufactured.T, strict=True)
            },
            "E_abs_opalx_V_per_m": e_opalx_abs,
            "B_abs_opalx_T": b_opalx_abs,
            "E_abs_manufactured_V_per_m": e_manufactured_abs,
            "B_abs_manufactured_T": b_manufactured_abs,
            "E_opalx_finite": e_finite,
            "B_opalx_finite": b_finite,
        }
    )
    table.to_csv(args.case_dir / "primary_field_diagnostic.csv", index=False)

    summary = {
        "green_function": green_function,
        "mesh": mesh,
        "probe_grid": [parameters.nx, parameters.nz],
        "probe_count": expected,
        "source_centroid_from_ip_m": source_center_from_ip_m,
        "manufactured_source_cutoff_sigma": 3.0,
        "finite_E_probe_fraction": float(np.mean(e_finite)),
        "finite_B_probe_fraction": float(np.mean(b_finite)),
        "max_manufactured_E_V_per_m": float(np.max(e_manufactured_abs)),
        "max_manufactured_B_T": float(np.max(b_manufactured_abs)),
        "relative_l2_E": (
            float(np.linalg.norm(e_opalx - e_manufactured) / np.linalg.norm(e_manufactured))
            if np.all(e_finite)
            else None
        ),
        "relative_l2_B": (
            float(np.linalg.norm(b_opalx - b_manufactured) / np.linalg.norm(b_manufactured))
            if np.all(b_finite)
            else None
        ),
    }
    (args.case_dir / "primary_field_diagnostic_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    cache = args.case_dir / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    (cache / "matplotlib").mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    shape = (parameters.nz, parameters.nx)
    z_grid = (positions_from_ip[:, 2] / parameters.sigma_z_m).reshape(shape)
    x_grid = (positions_from_ip[:, 0] / parameters.sigma_x_m).reshape(shape)
    quantities = (
        ("E", e_manufactured_abs, e_opalx_abs, e_finite, "V/m"),
        ("B", b_manufactured_abs, b_opalx_abs, b_finite, "T"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.8), sharex=True, sharey=True, constrained_layout=True)
    for row, (symbol, manufactured, opalx, finite, unit) in enumerate(quantities):
        vmin, vmax = positive_limits(manufactured)
        image = axes[row, 0].pcolormesh(
            z_grid,
            x_grid,
            manufactured.reshape(shape),
            shading="auto",
            cmap="magma",
            norm=LogNorm(vmin=vmin, vmax=vmax),
            rasterized=True,
        )
        axes[row, 0].set_title(rf"Manufactured $|{symbol}|$ [{unit}]")
        fig.colorbar(image, ax=axes[row, 0], pad=0.02)

        opalx_grid = np.ma.masked_invalid(opalx.reshape(shape))
        if np.any(finite):
            image = axes[row, 1].pcolormesh(
                z_grid,
                x_grid,
                opalx_grid,
                shading="auto",
                cmap="magma",
                norm=LogNorm(vmin=vmin, vmax=vmax),
                rasterized=True,
            )
            fig.colorbar(image, ax=axes[row, 1], pad=0.02)
        else:
            axes[row, 1].set_facecolor("0.85")
            axes[row, 1].text(0.5, 0.5, "no finite samples", ha="center", va="center", transform=axes[row, 1].transAxes)
        axes[row, 1].set_title(rf"OPALX $|{symbol}|$ [{unit}]")

        image = axes[row, 2].pcolormesh(
            z_grid,
            x_grid,
            finite.reshape(shape).astype(float),
            shading="auto",
            cmap="RdYlGn",
            vmin=0.0,
            vmax=1.0,
            rasterized=True,
        )
        axes[row, 2].set_title(rf"OPALX ${symbol}$ finite mask ({np.mean(finite):.1%})")
        fig.colorbar(image, ax=axes[row, 2], pad=0.02, ticks=(0.0, 1.0))
        axes[row, 0].set_ylabel(r"$x/\sigma_x$")

    for axis in axes[-1, :]:
        axis.set_xlabel(r"$(z-z_\mathrm{IP})/\sigma_z$")
    fig.suptitle(
        f"Primary-only BeamBeam field diagnostic: "
        f"{mesh[0]} x {mesh[1]} x {mesh[2]} {green_function.lower()} Green function; "
        r"manufactured source truncated at $3\sigma$"
    )
    for suffix in ("png", "pdf"):
        fig.savefig(args.case_dir / f"primary_field_diagnostic.{suffix}", dpi=240)
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
