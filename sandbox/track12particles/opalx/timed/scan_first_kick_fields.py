#!/usr/bin/env python3
"""Compare pair-4 OPALX probe fields with a rigid Gaussian reference.

The scan holds the exact-overlap track12 pair-4 geometry fixed and varies the
transverse field domain and mesh. It runs both the normal seeded OPALX Gaussian
sampler and a deterministic equal-probability tensor cubature. The latter is a
low-noise deposition diagnostic, not a replacement production distribution.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
TEMPLATE = HERE / "track12_pair4_field_probe.in.template"
DEFAULT_CONFIG = HERE / "first_kick_field_scan.json"
DEFAULT_OUTPUT = HERE / "first_kick_field_scan"
DEFAULT_OPALX = ROOT / "build_openmp" / "src" / "opalx"
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
REFERENCE = ROOT / "sandbox" / "TestParticleOrbit.dat"


def load_track12():
    spec = importlib.util.spec_from_file_location("track12_first_kick_model", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truncated_normal_midpoints(count: int, cutoff_sigma: float) -> np.ndarray:
    """Return symmetric midpoint quadrature nodes in truncated-normal probability."""
    if count <= 0 or cutoff_sigma <= 0.0:
        raise ValueError("tensor counts and cutoff must be positive")
    normal = statistics.NormalDist()
    lower = normal.cdf(-cutoff_sigma)
    upper = normal.cdf(+cutoff_sigma)
    probabilities = lower + (np.arange(count) + 0.5) * (upper - lower) / count
    values = np.array([normal.inv_cdf(float(value)) for value in probabilities])
    values -= np.mean(values)
    return values


def write_tensor_primary(
    path: Path, shape: tuple[int, int, int], cutoff_sigma: float, track12
) -> dict[str, object]:
    axes = [truncated_normal_midpoints(count, cutoff_sigma) for count in shape]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    positions = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    positions *= np.array([track12.SIGMA_XY_M, track12.SIGMA_XY_M, track12.SIGMA_Z_M])
    positions -= np.mean(positions, axis=0, dtype=np.float64)

    gamma = 1.0 + track12.DEFAULT_SOURCE_KINETIC_MEV * 1.0e6 / track12.ELECTRON_REST_EV
    beta_gamma = math.sqrt(gamma * gamma - 1.0)
    momenta = np.zeros_like(positions)
    momenta[:, 2] = beta_gamma

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(positions)}\n")
        stream.write("x y z px py pz\n")
        np.savetxt(stream, np.column_stack((positions, momenta)), fmt="%.16e")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "method": "equal-probability midpoint tensor cubature",
        "shape": list(shape),
        "particle_count": int(len(positions)),
        "cutoff_sigma": cutoff_sigma,
        "centroid_m": np.mean(positions, axis=0).tolist(),
        "beta_gamma": beta_gamma,
    }


def write_probe(path: Path, track12) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write("1\n")
        stream.write("x y z px py pz\n")
        stream.write(
            f"{track12.SIGMA_XY_M:.16e} 0.0000000000000000e+00 "
            "0.0000000000000000e+00 0.0000000000000000e+00 "
            "0.0000000000000000e+00 1.7320508075688772e+00\n"
        )


def primary_distribution(model: str) -> tuple[str, str]:
    if model == "random":
        return (
            """DIST_PRIMARY_ELECTRONS: DISTRIBUTION,
    TYPE = GAUSS,
    SIGMAX = primary_sigma_xy,
    SIGMAY = primary_sigma_xy,
    SIGMAZ = primary_sigma_z,
    SIGMAPX = 1.0e-12,
    SIGMAPY = 1.0e-12,
    SIGMAPZ = 1.0e-12,
    CUTOFFX = 6.0,
    CUTOFFY = 6.0,
    CUTOFFLONG = 6.0,
    NPARTDIST = primary_macroparticles;""",
            "    PC = primary_p0,",
        )
    if model == "tensor":
        return (
            """DIST_PRIMARY_ELECTRONS: DISTRIBUTION,
    TYPE = FROMFILE,
    FNAME = \"../input/primary_tensor.fromfile\",
    NPARTDIST = primary_macroparticles;""",
            "",
        )
    if model == "fixed":
        return (
            """DIST_PRIMARY_ELECTRONS: DISTRIBUTION,
    TYPE = FROMFILE,
    FNAME = \"../input/primary_fixed.fromfile\",
    NPARTDIST = primary_macroparticles;""",
            "",
        )
    raise ValueError(f"unsupported primary model {model!r}")


def render_case(
    case: dict[str, object], output_dir: Path, config: dict, tensor_count: int,
    fixed_primary: Path, track12
) -> Path:
    case_dir = output_dir / str(case["name"])
    input_dir = case_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    write_probe(input_dir / "pair4_probe.fromfile", track12)

    model = str(case["primary_model"])
    if model == "random":
        particle_count = int(config["primary"]["random_macroparticles"])
    elif model == "tensor":
        particle_count = tensor_count
    elif model == "fixed":
        if not fixed_primary.is_file():
            raise FileNotFoundError(fixed_primary)
        with fixed_primary.open("r", encoding="utf-8") as stream:
            particle_count = int(stream.readline().strip())
    else:
        raise ValueError(f"unsupported primary model {model!r}")
    distribution, beam_energy = primary_distribution(model)
    replacements = {
        "@PRIMARY_MACROPARTICLES@": str(particle_count),
        "@APERTURE_SPEC@": str(case["aperture"]),
        "@NX@": str(int(case["nx"])),
        "@NY@": str(int(case["ny"])),
        "@NZ@": str(int(config["mesh_nz"])),
        "@PRIMARY_DISTRIBUTION@": distribution,
        "@PRIMARY_BEAM_ENERGY@": beam_energy,
    }
    deck = TEMPLATE.read_text(encoding="utf-8")
    for token, value in replacements.items():
        deck = deck.replace(token, value)
    unresolved = [token for token in replacements if token in deck]
    if unresolved:
        raise RuntimeError(f"unresolved input tokens: {unresolved}")
    input_path = case_dir / "pair4_field_probe.in"
    input_path.write_text(deck, encoding="utf-8")
    return input_path


def run_case(opalx: Path, input_path: Path, omp_threads: int, force: bool) -> tuple[str, float]:
    h5_path = input_path.with_name(input_path.stem + "_c1.h5")
    if h5_path.exists() and not force:
        return "cached", 0.0
    start = time.monotonic()
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(omp_threads)
    completed = subprocess.run(
        [str(opalx), input_path.name],
        cwd=input_path.parent,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    runtime = time.monotonic() - start
    input_path.with_suffix(".log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"OPALX failed for {input_path.parent.name} with code {completed.returncode}; "
            f"see {input_path.with_suffix('.log')}"
        )
    return "ok", runtime


def external_runtime(input_path: Path) -> float:
    """Read the `/usr/bin/time -p` real time recorded by the A100 launcher."""
    path = input_path.parent / "runtime_compute.txt"
    if not path.is_file():
        return 0.0
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0] == "real":
            return float(fields[1])
    raise ValueError(f"{path}: missing real runtime")


def sorted_steps(h5_file: h5py.File) -> list[str]:
    return sorted(h5_file.keys(), key=lambda name: int(name.split("#", 1)[1]))


def read_probe(input_path: Path) -> dict[str, float]:
    h5_path = input_path.with_name(input_path.stem + "_c1.h5")
    with h5py.File(h5_path, "r") as h5_file:
        steps = sorted_steps(h5_file)
        if len(steps) < 2:
            raise RuntimeError(f"{h5_path}: expected at least two steps, found {steps}")
        group = h5_file[steps[1]]
        required = ("x", "y", "z", "Ex", "Ey", "Ez", "Bx", "By", "Bz")
        missing = [name for name in required if name not in group]
        if missing:
            raise KeyError(f"{h5_path}:{steps[1]} missing {missing}")
        if len(group["x"]) != 1:
            raise RuntimeError(f"{h5_path}:{steps[1]} expected one probe")
        ref_r = np.asarray(group.attrs.get("RefPartR", [0.0, 0.0, 0.0]), dtype=float)
        return {
            "x_m": float(group["x"][0] + ref_r[0]),
            "y_m": float(group["y"][0] + ref_r[1]),
            "z_abs_m": float(group["z"][0] + ref_r[2]),
            "Ex_opalx_V_per_m": float(group["Ex"][0]),
            "Ey_opalx_V_per_m": float(group["Ey"][0]),
            "Ez_opalx_V_per_m": float(group["Ez"][0]),
            "Bx_opalx_T": float(group["Bx"][0]),
            "By_opalx_T": float(group["By"][0]),
            "Bz_opalx_T": float(group["Bz"][0]),
        }


def source_centroid(input_path: Path) -> float:
    """Read absolute source z centroid from the second primary STAT row."""
    stat_paths = list(input_path.parent.glob("*_c0.stat"))
    if len(stat_paths) != 1:
        raise RuntimeError(f"{input_path.parent}: expected one primary stat file")
    text = stat_paths[0].read_text(encoding="utf-8").splitlines()
    column_names: list[str] = []
    parameter_count = 0
    data_start = None
    block = None
    for index, raw in enumerate(text):
        line = raw.strip()
        if line.startswith("&parameter"):
            parameter_count += 1
            block = "parameter"
        elif line.startswith("&column"):
            block = "column"
        elif line.startswith("&data"):
            block = "data"
        elif line.startswith("&end"):
            if block == "data":
                data_start = index + 1
                break
            block = None
        elif block == "column":
            match = re.match(r"name\s*=\s*([^,\s]+)", line)
            if match:
                column_names.append(match.group(1))
    if data_start is None:
        raise RuntimeError(f"{stat_paths[0]}: missing SDDS data block")
    rows = [line.split() for line in text[data_start:] if line.strip()][parameter_count:]
    if len(rows) < 2:
        raise RuntimeError(f"{stat_paths[0]}: expected at least two data rows")
    row = dict(zip(column_names, map(float, rows[1]), strict=True))
    return row["ref_z"] + row["mean_s"]


def manufactured_fields(
    position: np.ndarray,
    source_z_from_ip: float,
    track12,
    cutoff_sigma: float | None,
):
    beta = track12.beta_from_kinetic_energy(track12.DEFAULT_SOURCE_KINETIC_MEV)
    kwargs = {
        "charge_C": track12.DEFAULT_SOURCE_CHARGE_C,
        "sigma_lab_m": (track12.SIGMA_XY_M, track12.SIGMA_XY_M, track12.SIGMA_Z_M),
        "t0_s": 0.0,
        "cutoff_sigma": cutoff_sigma,
    }
    physical = track12.RigidAnisotropicGaussianSource(
        "physical", beta_z=+beta, center_t0_m=np.array([0.0, 0.0, source_z_from_ip]), **kwargs
    )
    copied = track12.RigidAnisotropicGaussianSource(
        "copied", beta_z=-beta, center_t0_m=np.array([0.0, 0.0, -source_z_from_ip]), **kwargs
    )
    return track12.anisotropic_total_lab_fields(position, 0.0, (physical, copied))


def interaction_point_s(input_path: Path) -> float:
    """Derive the IP from the placed BeamBeam element midpoint."""
    deck = input_path.read_text(encoding="utf-8")
    length_match = re.search(r"REAL\s+bb_length\s*=\s*([^;]+);", deck)
    edge_match = re.search(r"IP1\s*:\s*BEAMBEAM,.*?ELEMEDGE\s*=\s*([^,;]+)", deck, re.S)
    if length_match is None or edge_match is None:
        raise ValueError(f"{input_path}: cannot derive BeamBeam midpoint")
    return float(edge_match.group(1)) + 0.5 * float(length_match.group(1))


def cain_pair4_first_kick(track12) -> float:
    reference = track12.parse_reference_file(REFERENCE)
    rows = reference.loc[reference["pair"].eq(4) & reference["kind"].eq(2)].sort_values("step")
    return float((rows.iloc[1]["Px"] - rows.iloc[0]["Px"]) / track12.ELECTRON_REST_EV)


def field_kick(e_field: np.ndarray, b_field: np.ndarray, track12) -> float:
    """Constant-field Boris kick over one 1.8 um/c CAIN output interval."""
    initial = np.array([0.0, 0.0, math.sqrt(3.0)])
    dt_s = 1.8e-6 / track12.C_LIGHT
    final = track12.boris_kick(initial, e_field, b_field, dt_s, charge_units=-1.0)
    return float(final[0] - initial[0])


def analyze_case(
    case: dict[str, object], input_path: Path, status: str, runtime_s: float, track12
) -> dict[str, object]:
    row: dict[str, object] = {
        "case": str(case["name"]),
        "primary_model": str(case["primary_model"]),
        "aperture": str(case["aperture"]),
        "full_x_m": float(case["full_x_m"]),
        "full_y_m": float(case["full_y_m"]),
        "nx": int(case["nx"]),
        "ny": int(case["ny"]),
        "status": status,
        "runtime_s": runtime_s,
    }
    row["cell_x_um"] = float(case["full_x_m"]) / (int(case["nx"]) - 1) * 1.0e6
    row["cell_y_um"] = float(case["full_y_m"]) / (int(case["ny"]) - 1) * 1.0e6
    row.update(read_probe(input_path))

    ip_s = interaction_point_s(input_path)
    source_z = source_centroid(input_path) - ip_s
    position = np.array([row["x_m"], row["y_m"], row["z_abs_m"] - ip_s], dtype=float)
    e_analytic, b_analytic = manufactured_fields(position, source_z, track12, None)
    cutoff_sigma = float(case.get("manufactured_cutoff_sigma", 3.0))
    e_truncated, b_truncated = manufactured_fields(
        position, source_z, track12, cutoff_sigma
    )
    e_opalx = np.array(
        [row["Ex_opalx_V_per_m"], row["Ey_opalx_V_per_m"], row["Ez_opalx_V_per_m"]]
    )
    b_opalx = np.array([row["Bx_opalx_T"], row["By_opalx_T"], row["Bz_opalx_T"]])
    for axis, value in zip("xyz", e_analytic):
        row[f"E{axis}_analytic_V_per_m"] = float(value)
    for axis, value in zip("xyz", b_analytic):
        row[f"B{axis}_analytic_T"] = float(value)
    for axis, value in zip("xyz", e_truncated):
        row[f"E{axis}_truncated_V_per_m"] = float(value)
    for axis, value in zip("xyz", b_truncated):
        row[f"B{axis}_truncated_T"] = float(value)
    row["manufactured_cutoff_sigma"] = cutoff_sigma
    row["E_abs_opalx_V_per_m"] = float(np.linalg.norm(e_opalx))
    row["E_abs_analytic_V_per_m"] = float(np.linalg.norm(e_analytic))
    row["B_abs_opalx_T"] = float(np.linalg.norm(b_opalx))
    row["B_abs_analytic_T"] = float(np.linalg.norm(b_analytic))
    row["E_abs_truncated_V_per_m"] = float(np.linalg.norm(e_truncated))
    row["B_abs_truncated_T"] = float(np.linalg.norm(b_truncated))
    row["E_magnitude_ratio"] = row["E_abs_opalx_V_per_m"] / row["E_abs_analytic_V_per_m"]
    row["E_direction_cosine"] = float(
        np.dot(e_opalx, e_analytic)
        / (np.linalg.norm(e_opalx) * np.linalg.norm(e_analytic))
    )
    row["E_magnitude_ratio_truncated"] = (
        row["E_abs_opalx_V_per_m"] / row["E_abs_truncated_V_per_m"]
    )
    row["E_direction_cosine_truncated"] = float(
        np.dot(e_opalx, e_truncated)
        / (np.linalg.norm(e_opalx) * np.linalg.norm(e_truncated))
    )
    row["dPx_field_opalx"] = field_kick(e_opalx, b_opalx, track12)
    row["dPx_field_analytic"] = field_kick(e_analytic, b_analytic, track12)
    row["dPx_field_truncated"] = field_kick(e_truncated, b_truncated, track12)
    row["dPx_CAIN"] = cain_pair4_first_kick(track12)
    row["dPx_opalx_over_analytic"] = row["dPx_field_opalx"] / row["dPx_field_analytic"]
    row["dPx_opalx_over_truncated"] = (
        row["dPx_field_opalx"] / row["dPx_field_truncated"]
    )
    row["dPx_opalx_over_CAIN"] = row["dPx_field_opalx"] / row["dPx_CAIN"]
    return row


def configure_matplotlib(output_dir: Path) -> None:
    cache = output_dir / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def fitted_error_scaling(
    cell_spacing_um: pd.Series, ratio: pd.Series, samples: int = 100
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit ``abs(1-ratio) = C h**p`` and return the ratio curve and order."""
    h = np.asarray(cell_spacing_um, dtype=float)
    values = np.asarray(ratio, dtype=float)
    error = np.abs(1.0 - values)
    valid = np.isfinite(h) & np.isfinite(error) & (h > 0.0) & (error > 0.0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("a scaling fit requires at least two finite nonzero errors")
    order, log_coefficient = np.polyfit(np.log(h[valid]), np.log(error[valid]), 1)
    h_fit = np.geomspace(np.min(h[valid]), np.max(h[valid]), samples)
    side = np.sign(np.median(values[valid] - 1.0))
    if side == 0.0:
        side = -1.0
    ratio_fit = 1.0 + side * np.exp(log_coefficient) * h_fit**order
    return h_fit, ratio_fit, float(order)


def plot_results(results: pd.DataFrame, output_dir: Path) -> None:
    configure_matplotlib(output_dir)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)
    markers = {"random": "o", "tensor": "s", "fixed": "^"}
    refinement = results.loc[
        results["case"].str.match(r"square20_n(?:32|64|128)_(?:random|tensor)$")
    ]
    plotted = refinement if not refinement.empty else results
    for model, group in plotted.groupby("primary_model", sort=False):
        group = group.sort_values("cell_x_um")
        plot = axes[0].plot if len(group) > 1 else axes[0].scatter
        kwargs = {
            "marker": markers[model],
            "label": f"{model} source",
        }
        if len(group) > 1:
            kwargs.update({"lw": 1.5, "ms": 6.5})
        plot(group["cell_x_um"], group["E_magnitude_ratio_truncated"], **kwargs)

        plot = axes[1].plot if len(group) > 1 else axes[1].scatter
        plot(group["cell_x_um"], group["dPx_opalx_over_truncated"], **kwargs)

        proportional = group.loc[
            np.isclose(group["nx"] / group["ny"], 8.0)
            & group["case"].str.startswith("production_rect_")
        ].sort_values("cell_x_um", ascending=False)
        if model == "fixed" and len(proportional) >= 7:
            fit_groups = (
                ("3 coarse points", proportional.iloc[:3], "--", "C1"),
                ("4 fine points", proportional.iloc[-4:], "-.", "C3"),
            )
            for axis, ratio_column in zip(
                axes,
                ("E_magnitude_ratio_truncated", "dPx_opalx_over_truncated"),
            ):
                for label, fit_group, linestyle, color in fit_groups:
                    h_fit, ratio_fit, order = fitted_error_scaling(
                        fit_group["cell_x_um"], fit_group[ratio_column]
                    )
                    axis.plot(
                        h_fit,
                        ratio_fit,
                        linestyle=linestyle,
                        color=color,
                        lw=1.5,
                        label=rf"{label}: $1-r\propto h^{{{order:.2f}}}$",
                    )
    for axis in axes:
        axis.axhline(1.0, color="0.2", ls="--", lw=1.0)
        axis.invert_xaxis()
        axis.grid(True, which="both", color="0.88", lw=0.6)
        axis.set_xlabel(r"$x$ cell spacing [$\mu$m]")
        axis.legend(title="data and scaling fits", fontsize=8)
    axes[0].set_ylabel(r"$|E|_{\rm OPALX}/|E|_{\rm truncated}$")
    axes[1].set_ylabel(r"$\Delta p_x^{\rm OPALX}/\Delta p_x^{\rm truncated}$")
    if not refinement.empty:
        fig.suptitle(r"Pair-4 overlap: $40\,\mu$m full square-domain refinement")
    elif plotted["case"].str.startswith("production_rect_").all():
        fig.suptitle(
            r"Pair-4 overlap: production aperture, fixed $N_z=128$, $\pm3\sigma$ source"
        )
    fig.savefig(output_dir / "pair4_field_and_kick_convergence.png", dpi=240)
    plt.close(fig)


def add_observed_orders(results: pd.DataFrame) -> pd.DataFrame:
    """Add pairwise convergence orders for fixed-domain refinement sequences."""
    output = results.copy()
    output["effective_cell_um"] = np.sqrt(
        output["cell_x_um"] * output["cell_y_um"]
    )
    output["observed_order_E"] = np.nan
    output["observed_order_kick"] = np.nan
    group_columns = ["primary_model", "full_x_m", "full_y_m"]
    for _key, group in output.groupby(group_columns, sort=False):
        indices = group.sort_values("effective_cell_um", ascending=False).index
        for coarse_index, fine_index in zip(indices[:-1], indices[1:]):
            h_ratio = (
                output.at[coarse_index, "effective_cell_um"]
                / output.at[fine_index, "effective_cell_um"]
            )
            for ratio_column, order_column in (
                ("E_magnitude_ratio_truncated", "observed_order_E"),
                ("dPx_opalx_over_truncated", "observed_order_kick"),
            ):
                coarse_error = abs(1.0 - output.at[coarse_index, ratio_column])
                fine_error = abs(1.0 - output.at[fine_index, ratio_column])
                if h_ratio > 1.0 and coarse_error > 0.0 and fine_error > 0.0:
                    output.at[fine_index, order_column] = math.log(
                        coarse_error / fine_error
                    ) / math.log(h_ratio)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument(
        "--fixed-primary",
        type=Path,
        default=HERE / "input" / "primary_fixed.fromfile",
    )
    parser.add_argument("--omp-threads", type=int, default=8)
    parser.add_argument("--case", action="append", help="Run only named case(s)")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--max-e-magnitude-relative-error",
        type=float,
        help="Fail if any selected case exceeds this |E|-ratio error from unity.",
    )
    parser.add_argument(
        "--min-e-direction-cosine",
        type=float,
        help="Fail if any selected case has a smaller E direction cosine.",
    )
    args = parser.parse_args()

    if args.prepare_only and args.analyze_only:
        raise ValueError("--prepare-only and --analyze-only are mutually exclusive")
    if not args.analyze_only and not args.prepare_only and not args.opalx.is_file():
        raise FileNotFoundError(args.opalx)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    cases = config["cases"]
    if args.case:
        wanted = set(args.case)
        cases = [case for case in cases if case["name"] in wanted]
        missing = wanted - {case["name"] for case in cases}
        if missing:
            raise ValueError(f"unknown cases: {sorted(missing)}")
    track12 = load_track12()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fixed_primary = args.fixed_primary.resolve()
    if any(case["primary_model"] == "fixed" for case in cases) and not args.analyze_only:
        if not fixed_primary.is_file():
            raise FileNotFoundError(fixed_primary)
        staged_fixed_primary = args.output_dir / "input" / "primary_fixed.fromfile"
        staged_fixed_primary.parent.mkdir(parents=True, exist_ok=True)
        if staged_fixed_primary.resolve() != fixed_primary:
            shutil.copy2(fixed_primary, staged_fixed_primary)
        fixed_primary = staged_fixed_primary.resolve()

    tensor_metadata = None
    if any(case["primary_model"] == "tensor" for case in cases):
        tensor_shape = tuple(int(value) for value in config["primary"]["tensor_shape"])
        tensor_metadata = write_tensor_primary(
            args.output_dir / "input" / "primary_tensor.fromfile",
            tensor_shape,
            float(config["primary"]["cutoff_sigma"]),
            track12,
        )
        (args.output_dir / "input" / "primary_tensor_metadata.json").write_text(
            json.dumps(tensor_metadata, indent=2) + "\n", encoding="utf-8"
        )

    rows = []
    for case in cases:
        if args.analyze_only:
            input_path = (
                args.output_dir / str(case["name"]) / "pair4_field_probe.in"
            )
            if not input_path.is_file():
                raise FileNotFoundError(input_path)
        else:
            input_path = render_case(
                case,
                args.output_dir,
                config,
                int(tensor_metadata["particle_count"]) if tensor_metadata else 0,
                fixed_primary,
                track12,
            )
        if args.prepare_only:
            print(f"prepared {input_path}")
            continue
        if args.analyze_only:
            status, runtime = "external", external_runtime(input_path)
        else:
            status, runtime = run_case(
                args.opalx.resolve(), input_path, args.omp_threads, args.force
            )
        row = analyze_case(case, input_path, status, runtime, track12)
        rows.append(row)
        print(
            f"{row['case']}: |E| ratio={row['E_magnitude_ratio_truncated']:.6g}, "
            f"kick ratio={row['dPx_opalx_over_truncated']:.6g}, runtime={runtime:.2f} s",
            flush=True,
        )

    if args.prepare_only:
        return 0

    results = add_observed_orders(pd.DataFrame(rows))
    results.to_csv(args.output_dir / "pair4_first_kick_field_scan.csv", index=False)
    if not args.no_plot:
        plot_results(results, args.output_dir)
    external_manifests = []
    if args.analyze_only:
        for case in cases:
            run_manifest = args.output_dir / str(case["name"]) / "run_manifest.txt"
            if run_manifest.is_file():
                external_manifests.append(
                    {
                        "path": str(run_manifest.resolve()),
                        "sha256": sha256(run_manifest),
                    }
                )
    manifest = {
        "config": str(args.config.resolve()),
        "execution": "external" if args.analyze_only else "local",
        "opalx": (
            {
                "path": str(args.opalx.resolve()),
                "sha256": sha256(args.opalx.resolve()),
            }
            if not args.analyze_only
            else None
        ),
        "external_run_manifests": external_manifests,
        "tensor_primary": tensor_metadata,
        "fixed_primary": (
            {
                "path": str(args.fixed_primary.resolve()),
                "sha256": sha256(args.fixed_primary.resolve()),
            }
            if any(case["primary_model"] == "fixed" for case in cases)
            else None
        ),
        "result_rows": len(results),
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    columns = [
        "case",
        "E_magnitude_ratio",
        "E_magnitude_ratio_truncated",
        "E_direction_cosine_truncated",
        "dPx_field_opalx",
        "dPx_field_analytic",
        "dPx_field_truncated",
        "dPx_CAIN",
        "observed_order_E",
        "observed_order_kick",
    ]
    print(results[columns].to_string(index=False))
    print(f"wrote {args.output_dir}")

    failures = []
    if args.max_e_magnitude_relative_error is not None:
        errors = np.abs(results["E_magnitude_ratio_truncated"] - 1.0)
        failed = results.loc[errors.gt(args.max_e_magnitude_relative_error), "case"]
        if not failed.empty:
            failures.append(
                "|E| relative error exceeded "
                f"{args.max_e_magnitude_relative_error:g}: {failed.tolist()}"
            )
    if args.min_e_direction_cosine is not None:
        failed = results.loc[
            results["E_direction_cosine_truncated"].lt(args.min_e_direction_cosine), "case"
        ]
        if not failed.empty:
            failures.append(
                "E direction cosine fell below "
                f"{args.min_e_direction_cosine:g}: {failed.tolist()}"
            )
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
