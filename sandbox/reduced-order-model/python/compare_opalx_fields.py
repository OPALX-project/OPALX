#!/usr/bin/env python3
"""Compare OPALX rigid-source probe fields with the smooth Gaussian model."""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent
FIELD_MODEL_SCRIPT = SCRIPT_DIR / "rigid_two_gaussian_fields.py"
DEFAULT_CASE_DIR = MODEL_DIR / "outputs" / "opalx-fields"
DEFAULT_OUTPUT_DIR = MODEL_DIR / "outputs" / "comparison"


def load_field_model():
    spec = importlib.util.spec_from_file_location("rigid_two_gaussian_fields", FIELD_MODEL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {FIELD_MODEL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sorted_steps(h5_file: h5py.File) -> list[str]:
    return sorted(h5_file.keys(), key=lambda key: int(key.split("#", maxsplit=1)[1]))


def reference_z(group: h5py.Group) -> float:
    if "RefPartR" in group.attrs:
        return float(group.attrs["RefPartR"][2])
    if "SPOS" in group.attrs:
        return float(group.attrs["SPOS"][0])
    raise KeyError("H5 step has neither RefPartR nor SPOS")


def read_source_centroid_from_stat(path: Path, step_index: int) -> float:
    """Read ref_z + mean_s from an ASCII SDDS STAT file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    column_names = []
    parameter_count = 0
    data_end = None
    block = None
    for index, raw_line in enumerate(lines):
        line = raw_line.strip()
        if line.startswith("&parameter"):
            parameter_count += 1
            block = "parameter"
        elif line.startswith("&column"):
            block = "column"
        elif line.startswith("&data"):
            block = "data"
        elif line.startswith("&end"):
            if block == "data":
                data_end = index
                break
            block = None
        elif block == "column":
            match = re.match(r"name\s*=\s*([^,\s]+)", line)
            if match:
                column_names.append(match.group(1))

    if data_end is None:
        raise ValueError(f"{path}: missing SDDS data section")
    required = ("mean_s", "ref_z")
    missing = [name for name in required if name not in column_names]
    if missing:
        raise KeyError(f"{path}: missing STAT columns {missing}")

    data_lines = [line.strip() for line in lines[data_end + 1 :] if line.strip()]
    rows = data_lines[parameter_count:]
    if step_index < 0 or step_index >= len(rows):
        raise IndexError(f"{path}: requested row {step_index}, found {len(rows)} rows")
    values = rows[step_index].split()
    if len(values) != len(column_names):
        raise ValueError(
            f"{path}: row {step_index} has {len(values)} values for "
            f"{len(column_names)} columns"
        )
    row = dict(zip(column_names, map(float, values)))
    return row["ref_z"] + row["mean_s"]


def read_case(case_dir: Path, step_index: int) -> tuple[pd.DataFrame, float]:
    witness_files = sorted(case_dir.glob("*_c1.h5"))
    source_files = sorted(case_dir.glob("*_c0.h5"))
    source_stat_files = sorted(case_dir.glob("*_c0.stat"))
    if len(witness_files) != 1:
        raise FileNotFoundError(
            f"{case_dir}: expected one _c1.h5, found {len(witness_files)}"
        )

    if len(source_stat_files) == 1:
        source_centroid_abs_z_m = read_source_centroid_from_stat(
            source_stat_files[0], step_index
        )
    elif len(source_files) == 1:
        with h5py.File(source_files[0], "r") as source_h5:
            step_name = f"Step#{step_index}"
            if step_name not in source_h5:
                raise KeyError(f"{source_files[0]} has steps {sorted_steps(source_h5)}")
            group = source_h5[step_name]
            source_centroid_abs_z_m = reference_z(group) + float(np.mean(group["z"][:]))
    else:
        raise FileNotFoundError(
            f"{case_dir}: expected one _c0.stat or one _c0.h5, found "
            f"{len(source_stat_files)} and {len(source_files)}"
        )

    with h5py.File(witness_files[0], "r") as witness_h5:
        step_name = f"Step#{step_index}"
        if step_name not in witness_h5:
            raise KeyError(f"{witness_files[0]} has steps {sorted_steps(witness_h5)}")
        group = witness_h5[step_name]
        required = ("x", "y", "z", "Ex", "Ey", "Ez", "Bx", "By", "Bz")
        missing = [name for name in required if name not in group]
        if missing:
            raise KeyError(f"{witness_files[0]}:{step_name} missing {missing}; set EBDUMP=TRUE")
        frame = pd.DataFrame({name: group[name][:] for name in required})
        frame["z_abs_m"] = reference_z(group) + frame["z"]
        if "id" in group:
            frame["id"] = group["id"][:]
    return frame, source_centroid_abs_z_m


def actual_sources(
    field_model,
    track12,
    parameters,
    source_center_from_ip_m: float,
    source_model: str,
):
    beta = track12.beta_from_kinetic_energy(parameters.kinetic_energy_MeV)
    charge_C = -parameters.electrons_per_bunch * track12.ELEMENTARY_CHARGE_C
    sigma_lab_m = (parameters.sigma_x_m, parameters.sigma_y_m, parameters.sigma_z_m)
    physical = track12.RigidAnisotropicGaussianSource(
        "physical",
        charge_C=charge_C,
        sigma_lab_m=sigma_lab_m,
        beta_z=+beta,
        center_t0_m=np.array([0.0, 0.0, source_center_from_ip_m]),
    )
    copied = track12.RigidAnisotropicGaussianSource(
        "copied",
        charge_C=charge_C,
        sigma_lab_m=sigma_lab_m,
        beta_z=-beta,
        center_t0_m=np.array([0.0, 0.0, -source_center_from_ip_m]),
    )
    if source_model == "physical":
        return (physical,)
    if source_model == "physical-and-copied":
        return physical, copied
    raise ValueError(f"unknown source model {source_model!r}")


def analyze_case(
    case_dir: Path,
    field_model,
    track12,
    parameters,
    step_index: int,
    source_model: str = "physical-and-copied",
):
    frame, source_centroid_abs_z_m = read_case(case_dir, step_index)
    ip_s_m = 0.008
    source_center_from_ip_m = source_centroid_abs_z_m - ip_s_m
    positions = frame[["x", "y"]].to_numpy(dtype=float)
    positions = np.column_stack((positions, frame["z_abs_m"].to_numpy() - ip_s_m))
    sources = actual_sources(
        field_model, track12, parameters, source_center_from_ip_m, source_model
    )
    e_analytic, b_analytic = field_model.total_lab_fields_batch(
        positions, sources, track12, parameters.quadrature_chunk_size
    )
    e_opalx = frame[["Ex", "Ey", "Ez"]].to_numpy(dtype=float)
    b_opalx = frame[["Bx", "By", "Bz"]].to_numpy(dtype=float)

    output = pd.DataFrame(
        {
            "case": case_dir.name,
            "x_m": positions[:, 0],
            "y_m": positions[:, 1],
            "z_minus_ip_m": positions[:, 2],
            "Ex_opalx_V_per_m": e_opalx[:, 0],
            "Ey_opalx_V_per_m": e_opalx[:, 1],
            "Ez_opalx_V_per_m": e_opalx[:, 2],
            "Bx_opalx_T": b_opalx[:, 0],
            "By_opalx_T": b_opalx[:, 1],
            "Bz_opalx_T": b_opalx[:, 2],
            "Ex_analytic_V_per_m": e_analytic[:, 0],
            "Ey_analytic_V_per_m": e_analytic[:, 1],
            "Ez_analytic_V_per_m": e_analytic[:, 2],
            "Bx_analytic_T": b_analytic[:, 0],
            "By_analytic_T": b_analytic[:, 1],
            "Bz_analytic_T": b_analytic[:, 2],
        }
    )
    e_opalx_abs = np.linalg.norm(e_opalx, axis=1)
    e_analytic_abs = np.linalg.norm(e_analytic, axis=1)
    b_opalx_abs = np.linalg.norm(b_opalx, axis=1)
    b_analytic_abs = np.linalg.norm(b_analytic, axis=1)
    output["E_abs_opalx_V_per_m"] = e_opalx_abs
    output["E_abs_analytic_V_per_m"] = e_analytic_abs
    output["B_abs_opalx_T"] = b_opalx_abs
    output["B_abs_analytic_T"] = b_analytic_abs

    e_reference_nonzero = e_analytic_abs > np.max(e_analytic_abs) * 1.0e-12
    b_reference_nonzero = b_analytic_abs > max(
        np.max(b_analytic_abs) * 1.0e-12, 1.0e-14
    )
    e_opalx_nonzero = e_opalx_abs > np.max(e_opalx_abs) * 1.0e-12
    b_opalx_nonzero = b_opalx_abs > max(np.max(b_opalx_abs) * 1.0e-12, 1.0e-14)
    e_direction_valid = e_reference_nonzero & e_opalx_nonzero
    b_direction_valid = b_reference_nonzero & b_opalx_nonzero
    uncovered = (e_reference_nonzero & ~e_opalx_nonzero) | (
        b_reference_nonzero & ~b_opalx_nonzero
    )
    output["E_magnitude_ratio_opalx_over_analytic"] = np.divide(
        e_opalx_abs,
        e_analytic_abs,
        out=np.full_like(e_opalx_abs, np.nan),
        where=e_analytic_abs > np.max(e_analytic_abs) * 1.0e-12,
    )
    output["B_magnitude_ratio_opalx_over_analytic"] = np.divide(
        b_opalx_abs,
        b_analytic_abs,
        out=np.full_like(b_opalx_abs, np.nan),
        where=b_analytic_abs > max(np.max(b_analytic_abs) * 1.0e-12, 1.0e-14),
    )
    e_direction = np.sum(
        e_opalx[e_direction_valid] * e_analytic[e_direction_valid], axis=1
    ) / (
        e_opalx_abs[e_direction_valid] * e_analytic_abs[e_direction_valid]
    )
    b_direction = (
        np.sum(b_opalx[b_direction_valid] * b_analytic[b_direction_valid], axis=1)
        / (b_opalx_abs[b_direction_valid] * b_analytic_abs[b_direction_valid])
        if np.any(b_direction_valid)
        else np.array([], dtype=float)
    )
    summary = {
        "case": case_dir.name,
        "source_model": source_model,
        "requested_separation_sigma_z": float(
            case_dir.name.replace("p", ".").replace("sigma", "")
        ),
        "actual_separation_sigma_z": 2.0 * abs(source_center_from_ip_m) / parameters.sigma_z_m,
        "field_probe_count": len(output),
        "max_E_opalx_V_per_m": float(np.max(e_opalx_abs)),
        "max_E_analytic_V_per_m": float(np.max(e_analytic_abs)),
        "max_B_opalx_T": float(np.max(b_opalx_abs)),
        "max_B_analytic_T": float(np.max(b_analytic_abs)),
        "relative_l2_E": float(np.linalg.norm(e_opalx - e_analytic) / np.linalg.norm(e_analytic)),
        "median_E_magnitude_ratio": float(
            np.median(e_opalx_abs[e_reference_nonzero] / e_analytic_abs[e_reference_nonzero])
        ),
        "median_E_direction_cosine": float(np.nanmedian(e_direction)),
        "near_zero_E_fraction_opalx": float(1.0 - np.mean(e_opalx_nonzero)),
        "relative_l2_B": (
            float(np.linalg.norm(b_opalx - b_analytic) / np.linalg.norm(b_analytic))
            if np.linalg.norm(b_analytic) > 0.0
            else float("nan")
        ),
        "median_B_magnitude_ratio": (
            float(
                np.median(
                    b_opalx_abs[b_reference_nonzero] / b_analytic_abs[b_reference_nonzero]
                )
            )
            if np.any(b_reference_nonzero)
            else float("nan")
        ),
        "median_B_direction_cosine": (
            float(np.nanmedian(b_direction)) if len(b_direction) else float("nan")
        ),
        "near_zero_B_fraction_opalx": float(1.0 - np.mean(b_opalx_nonzero)),
        "uncovered_probe_fraction": float(np.mean(uncovered)),
    }
    return output, summary


def configure_matplotlib(output_dir: Path) -> None:
    cache_root = output_dir / ".plot-cache"
    cache_dir = cache_root / "matplotlib"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    cache_dir.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def plot_case(
    comparison: pd.DataFrame,
    summary: dict[str, float],
    parameters,
    output_dir: Path,
) -> Path:
    configure_matplotlib(output_dir)
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.colors import LogNorm

    z_normalized = comparison["z_minus_ip_m"].to_numpy() / parameters.sigma_z_m
    x_normalized = comparison["x_m"].to_numpy() / parameters.sigma_x_m
    triangulation = mtri.Triangulation(z_normalized, x_normalized)
    quantities = (
        (
            "E",
            comparison["E_abs_analytic_V_per_m"].to_numpy(),
            comparison["E_abs_opalx_V_per_m"].to_numpy(),
            comparison["E_magnitude_ratio_opalx_over_analytic"].to_numpy(),
            "V/m",
        ),
        (
            "B",
            comparison["B_abs_analytic_T"].to_numpy(),
            comparison["B_abs_opalx_T"].to_numpy(),
            comparison["B_magnitude_ratio_opalx_over_analytic"].to_numpy(),
            "T",
        ),
    )
    fig, axes = plt.subplots(
        2, 3, figsize=(11.5, 6.9), sharex=True, sharey=True, constrained_layout=True
    )
    for row, (symbol, analytic, opalx, ratio, unit) in enumerate(quantities):
        for column, (values, title) in enumerate(
            ((analytic, "manufactured"), (opalx, "OPALX"))
        ):
            positive = values[values > 0.0]
            maximum = float(np.max(positive))
            minimum = max(float(np.quantile(positive, 0.01)), maximum * 1.0e-6)
            image = axes[row, column].tripcolor(
                triangulation,
                values,
                shading="gouraud",
                cmap="magma",
                norm=LogNorm(vmin=minimum, vmax=maximum),
                rasterized=True,
            )
            axes[row, column].set_title(f"$|{symbol}|$ {title} [{unit}]")
            fig.colorbar(image, ax=axes[row, column], pad=0.02)

        finite_ratio = ratio[np.isfinite(ratio) & (ratio > 0.0)]
        log_ratio = np.full_like(ratio, np.nan)
        valid = np.isfinite(ratio) & (ratio > 0.0)
        log_ratio[valid] = np.log10(ratio[valid])
        ratio_min, ratio_max = np.quantile(np.log10(finite_ratio), (0.01, 0.99))
        ratio_image = axes[row, 2].tripcolor(
            triangulation,
            log_ratio,
            shading="gouraud",
            cmap="viridis",
            vmin=float(ratio_min),
            vmax=float(ratio_max),
            rasterized=True,
        )
        axes[row, 2].set_title(
            rf"$\log_{{10}}(|{symbol}|_\mathrm{{OPALX}}/|{symbol}|_\mathrm{{manuf.}})$"
        )
        fig.colorbar(ratio_image, ax=axes[row, 2], pad=0.02)
        axes[row, 0].set_ylabel(r"$x/\sigma_x$")

    half_separation = 0.5 * summary["actual_separation_sigma_z"]
    for axis in axes.ravel():
        axis.plot(
            [-half_separation, +half_separation],
            [0.0, 0.0],
            marker="x",
            linestyle="none",
            color="cyan",
            markersize=5,
            markeredgewidth=1.0,
        )
        axis.axvline(0.0, color="white", linewidth=0.5, alpha=0.55)
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$(z-z_\mathrm{IP})/\sigma_z$")

    fig.suptitle(
        f"OPALX {summary['source_model']} field versus manufactured solution\n"
        f"d={summary['actual_separation_sigma_z']:.6f} sigma_z; "
        f"relative L2(E)={summary['relative_l2_E']:.3g}, "
        f"relative L2(B)={summary['relative_l2_B']:.3g}",
        fontsize=11,
    )
    source_label = str(summary["source_model"]).replace("-", "_")
    path = output_dir / f"{summary['case']}_{source_label}_field_comparison.png"
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument(
        "--source-model",
        choices=("physical", "physical-and-copied"),
        default="physical-and-copied",
        help="Manufactured sources to compare against for the selected OPALX step.",
    )
    parser.add_argument(
        "--cases", nargs="+", default=None, help="Optional case directory names, e.g. 3sigma."
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    field_model = load_field_model()
    parameters = field_model.load_parameters(args.config or field_model.DEFAULT_CONFIG)
    track12 = field_model.load_track12_module()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    comparisons = []
    summaries = []
    for case_dir in sorted(path for path in args.case_dir.iterdir() if path.is_dir()):
        if args.cases is not None and case_dir.name not in args.cases:
            continue
        if not list(case_dir.glob("*_c1.h5")):
            continue
        comparison, summary = analyze_case(
            case_dir, field_model, track12, parameters, args.step, args.source_model
        )
        comparisons.append(comparison)
        summaries.append(summary)
        comparison.to_csv(args.output_dir / f"{case_dir.name}_field_comparison.csv", index=False)
        plot_case(comparison, summary, parameters, args.output_dir)

    if not comparisons:
        raise FileNotFoundError(f"no completed OPALX cases below {args.case_dir}")
    combined = pd.concat(comparisons, ignore_index=True)
    combined.to_csv(args.output_dir / "opalx_vs_rigid_gaussian_fields.csv", index=False)
    summary_frame = pd.DataFrame(summaries).sort_values(
        "requested_separation_sigma_z", ascending=False
    )
    summary_frame.to_csv(args.output_dir / "opalx_vs_rigid_gaussian_summary.csv", index=False)
    print(summary_frame.to_string(index=False, float_format=lambda value: f"{value:.6e}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
