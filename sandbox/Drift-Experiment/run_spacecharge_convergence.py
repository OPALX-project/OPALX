#!/usr/bin/env python3
"""Run the drift space-charge convergence matrix."""

from __future__ import annotations

import argparse
import importlib.util
import math
import os
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "spacecharge_drift_30cm.in"
DEFAULT_OPALX = REPO / "build_openmp" / "src" / "opalx"
COMPARATOR_MODULE = REPO / "sandbox" / "compare-e-fields" / "compare_gaussian_pic_fields.py"

N_GRID_VALUES = (16, 32, 64)
NPPG_VALUES = (5, 8, 14)
SEED_VALUES = (42, 43, 44)
SIGMA_M = 1.0e-3
TOTAL_CHARGE_C = -1.0e-9
EPS0 = 8.8541878128e-12


def load_comparator_module():
    if not COMPARATOR_MODULE.exists():
        return LocalComparator
    spec = importlib.util.spec_from_file_location("compare_gaussian_pic_fields", COMPARATOR_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {COMPARATOR_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class LocalComparator:
    @staticmethod
    def vector_abs(frame: pd.DataFrame, prefix: str) -> np.ndarray:
        return np.sqrt(
            frame[f"{prefix}Ex"].to_numpy(float) ** 2
            + frame[f"{prefix}Ey"].to_numpy(float) ** 2
            + frame[f"{prefix}Ez"].to_numpy(float) ** 2
        )

    @staticmethod
    def gaussian_e_field(
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        charge: float,
        sigma: float,
        center: tuple[float, float, float],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        dx = np.asarray(x, dtype=float) - center[0]
        dy = np.asarray(y, dtype=float) - center[1]
        dz = np.asarray(z, dtype=float) - center[2]
        r2 = dx * dx + dy * dy + dz * dz
        r = np.sqrt(r2)
        field_prefactor = np.zeros_like(r)
        nonzero = r > 0.0
        scaled = r[nonzero] / (math.sqrt(2.0) * sigma)
        enclosed_factor = (
            np.vectorize(math.erf)(scaled)
            - math.sqrt(2.0 / math.pi) * (r[nonzero] / sigma) * np.exp(-0.5 * (r[nonzero] / sigma) ** 2)
        )
        field_prefactor[nonzero] = charge / (4.0 * math.pi * EPS0) * enclosed_factor / (r[nonzero] ** 3)
        return field_prefactor * dx, field_prefactor * dy, field_prefactor * dz

    @staticmethod
    def add_error_columns(
        comparison: pd.DataFrame,
        rel_floor: float,
        center: tuple[float, float, float],
    ) -> pd.DataFrame:
        result = comparison.copy()
        result["dEx"] = result["pic_Ex"] - result["ana_Ex"]
        result["dEy"] = result["pic_Ey"] - result["ana_Ey"]
        result["dEz"] = result["pic_Ez"] - result["ana_Ez"]
        result["dE_vec_abs"] = np.sqrt(result["dEx"] ** 2 + result["dEy"] ** 2 + result["dEz"] ** 2)
        result["dE_abs"] = result["pic_E_abs"] - result["ana_E_abs"]
        result["abs_dE_abs"] = result["dE_abs"].abs()
        result["rel_vec_err"] = result["dE_vec_abs"] / np.maximum(result["ana_E_abs"], rel_floor)
        result["rel_mag_err"] = result["dE_abs"] / np.maximum(result["ana_E_abs"], rel_floor)
        dx = result["x_m"] - center[0]
        dy = result["y_m"] - center[1]
        dz = result["z_m"] - center[2]
        result["r_m"] = np.sqrt(dx * dx + dy * dy + dz * dz)
        return result

    @staticmethod
    def compute_metrics(comparison: pd.DataFrame, rel_floor: float) -> pd.DataFrame:
        del rel_floor
        d2 = comparison["dEx"] ** 2 + comparison["dEy"] ** 2 + comparison["dEz"] ** 2
        pic2 = comparison["pic_Ex"] ** 2 + comparison["pic_Ey"] ** 2 + comparison["pic_Ez"] ** 2
        ana2 = comparison["ana_Ex"] ** 2 + comparison["ana_Ey"] ** 2 + comparison["ana_Ez"] ** 2
        mag_d2 = comparison["dE_abs"] ** 2

        metrics = {
            "n_particles": len(comparison),
            "vector_error_L1_sum": float(comparison["dE_vec_abs"].sum()),
            "vector_error_L2_global": float(np.sqrt(d2.sum())),
            "vector_error_Linf_max": float(comparison["dE_vec_abs"].max()),
            "vector_error_RMSE": float(np.sqrt(d2.mean())),
            "relative_vector_L2_vs_analytic": float(np.sqrt(d2.sum()) / np.sqrt(ana2.sum())),
            "relative_vector_L2_vs_pic": float(np.sqrt(d2.sum()) / np.sqrt(pic2.sum())),
            "mean_pointwise_relative_vector_error": float(comparison["rel_vec_err"].mean()),
            "median_pointwise_relative_vector_error": float(comparison["rel_vec_err"].median()),
            "p95_pointwise_relative_vector_error": float(comparison["rel_vec_err"].quantile(0.95)),
            "magnitude_error_RMSE": float(np.sqrt(mag_d2.mean())),
            "magnitude_error_Linf_max": float(comparison["abs_dE_abs"].max()),
            "relative_magnitude_L2_vs_analytic": float(np.sqrt(mag_d2.sum()) / np.sqrt((comparison["ana_E_abs"] ** 2).sum())),
            "mean_abs_pic_E": float(comparison["pic_E_abs"].mean()),
            "mean_abs_analytic_E": float(comparison["ana_E_abs"].mean()),
            "max_abs_pic_E": float(comparison["pic_E_abs"].max()),
            "max_abs_analytic_E": float(comparison["ana_E_abs"].max()),
        }
        for comp in ("x", "y", "z"):
            err = comparison[f"dE{comp}"]
            ana = comparison[f"ana_E{comp}"]
            metrics[f"E{comp}_error_mean"] = float(err.mean())
            metrics[f"E{comp}_error_RMSE"] = float(np.sqrt((err**2).mean()))
            metrics[f"E{comp}_error_Linf_max"] = float(err.abs().max())
            metrics[f"E{comp}_relative_L2_vs_analytic"] = float(np.sqrt((err**2).sum()) / np.sqrt((ana**2).sum()))
        return pd.DataFrame([metrics])


def parse_csv_ints(value: str) -> tuple[int, ...]:
    return tuple(int(item.strip()) for item in value.split(",") if item.strip())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx", type=Path, default=DEFAULT_OPALX)
    parser.add_argument("--output-dir", type=Path, default=HERE / "convergence_results")
    parser.add_argument("--n-grid-values", default=",".join(map(str, N_GRID_VALUES)))
    parser.add_argument("--nppg-values", default=",".join(map(str, NPPG_VALUES)))
    parser.add_argument("--seeds", default=",".join(map(str, SEED_VALUES)))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Replace existing output directory")
    parser.add_argument("--keep-raw", action="store_true", help="Keep large OPALX H5 particle dumps")
    parser.add_argument("--resume", action="store_true", help="Skip cases that already have metrics.csv")
    parser.add_argument("--sample-particles", type=int, default=50000)
    return parser.parse_args()


def replace_assignment(text: str, name: str, value: int) -> str:
    prefix = f"REAL {name} = "
    lines = []
    replaced = False
    for line in text.splitlines():
        if line.strip().startswith(prefix):
            lines.append(f"REAL {name} = {value};")
            replaced = True
        else:
            lines.append(line)
    if not replaced:
        raise RuntimeError(f"Could not find assignment for {name} in {TEMPLATE}")
    return "\n".join(lines) + "\n"


def replace_option_seed(text: str, seed: int) -> str:
    lines = []
    replaced = False
    for line in text.splitlines():
        if line.strip().startswith("OPTION, SEED"):
            lines.append(f"OPTION, SEED = {seed};")
            replaced = True
        elif line.strip().startswith("OPTION, PSDUMPFREQ"):
            lines.append("OPTION, PSDUMPFREQ = 0;")
        else:
            lines.append(line)
    if not replaced:
        raise RuntimeError(f"Could not find OPTION, SEED in {TEMPLATE}")
    return "\n".join(lines) + "\n"


def render_input(n_grid: int, nppg: int, seed: int, path: Path) -> int:
    text = TEMPLATE.read_text(encoding="utf-8")
    text = replace_assignment(text, "N_GRID", n_grid)
    text = replace_assignment(text, "NPPG", nppg)
    text = replace_option_seed(text, seed)
    path.write_text(text, encoding="utf-8")
    return n_grid * n_grid * n_grid * nppg


def run_opalx(opalx: Path, deck: Path, threads: int, log: Path) -> None:
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OPALX_SC_FIELD_H5_STEPS"] = "0"
    with log.open("w", encoding="utf-8") as out:
        subprocess.run(
            [str(opalx.resolve()), deck.name],
            cwd=deck.parent,
            env=env,
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
        )


def read_h5_attr_scalar(attrs, name: str, default: int = 0) -> int:
    value = attrs.get(name)
    if value is None:
        return default
    array = np.asarray(value)
    if array.size == 0:
        return default
    return int(array.flat[0])


def read_particle_field_h5(case_dir: Path) -> pd.DataFrame:
    h5_files = sorted(case_dir.glob("*_c0.h5")) or sorted(case_dir.glob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(f"No OPALX H5 particle dump found in {case_dir}")
    h5_path = h5_files[0]

    with h5py.File(h5_path, "r") as h5:
        def step_sort_key(name: str) -> int:
            if name.startswith("Step#"):
                try:
                    return int(name.split("#", 1)[1])
                except ValueError:
                    pass
            return 10**12

        step_names = [
            name
            for name, obj in h5.items()
            if isinstance(obj, h5py.Group)
            and all(field in obj for field in ("x", "y", "z", "Ex", "Ey", "Ez"))
        ]
        if not step_names:
            raise RuntimeError(f"No H5 step with x/y/z/Ex/Ey/Ez found in {h5_path}")
        step_name = min(step_names, key=step_sort_key)
        step = h5[step_name]
        global_track_step = read_h5_attr_scalar(step.attrs, "GlobalTrackStep", default=0)
        if global_track_step != 0:
            raise RuntimeError(
                f"Expected the space-charge H5 diagnostic at GlobalTrackStep=0 in {h5_path}, "
                f"but selected {step_name} with GlobalTrackStep={global_track_step}. "
                "Check OPALX_SC_FIELD_H5_STEPS and phase-space dump settings."
            )

        n_particles = len(step["x"])
        rank = step["rank"][:] if "rank" in step else np.zeros(n_particles, dtype=int)
        return pd.DataFrame(
            {
                "step": np.zeros(n_particles, dtype=int),
                "container": np.zeros(n_particles, dtype=int),
                "rank": rank,
                "local_index": np.arange(n_particles, dtype=int),
                "x_m": step["x"][:],
                "y_m": step["y"][:],
                "z_m": step["z"][:],
                "Ex_V_per_m": step["Ex"][:],
                "Ey_V_per_m": step["Ey"][:],
                "Ez_V_per_m": step["Ez"][:],
            }
        )


def save_case_radial_plots(
    comparison: pd.DataFrame,
    center: tuple[float, float, float],
    case_dir: Path,
    max_plot_points: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    coords = comparison[["x_m", "y_m", "z_m"]].to_numpy(float)
    fields_pic = comparison[["pic_Ex", "pic_Ey", "pic_Ez"]].to_numpy(float)
    fields_ana = comparison[["ana_Ex", "ana_Ey", "ana_Ez"]].to_numpy(float)
    rvec = coords - np.asarray(center, dtype=float)
    radius = np.linalg.norm(rvec, axis=1)
    rhat = np.zeros_like(rvec)
    nonzero = radius > 0.0
    rhat[nonzero] = rvec[nonzero] / radius[nonzero, None]

    pic_abs = np.linalg.norm(fields_pic, axis=1)
    ana_abs = np.linalg.norm(fields_ana, axis=1)
    pic_er = np.einsum("ij,ij->i", fields_pic, rhat)
    ana_er = np.einsum("ij,ij->i", fields_ana, rhat)
    rel_mag = (pic_abs - ana_abs) / np.maximum(ana_abs, 1.0e-30)

    rmax = float(np.quantile(radius, 0.995))
    edges = np.linspace(0.0, rmax, 35)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        selected = (radius >= lo) & (radius < hi)
        if int(selected.sum()) < 8:
            continue
        rows.append(
            {
                "r_mean_mm": 1.0e3 * float(radius[selected].mean()),
                "n": int(selected.sum()),
                "pic_abs_mean_V_per_m": float(pic_abs[selected].mean()),
                "ana_abs_mean_V_per_m": float(ana_abs[selected].mean()),
                "pic_er_mean_V_per_m": float(pic_er[selected].mean()),
                "ana_er_mean_V_per_m": float(ana_er[selected].mean()),
                "rel_mag_mean": float(rel_mag[selected].mean()),
                "rel_mag_median": float(np.median(rel_mag[selected])),
            }
        )
    profile = pd.DataFrame(rows)
    profile.to_csv(case_dir / "radial_profile.csv", index=False)

    if max_plot_points > 0 and len(comparison) > max_plot_points:
        rng = np.random.default_rng(20260703)
        plot_index = np.sort(rng.choice(len(comparison), size=max_plot_points, replace=False))
    else:
        plot_index = np.arange(len(comparison))
    radius_mm = 1.0e3 * radius[plot_index]

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "figure.dpi": 140,
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ax.scatter(
        radius_mm,
        pic_abs[plot_index],
        s=7,
        alpha=0.16,
        color="#386cb0",
        linewidths=0,
        label="OPALX particles",
    )
    ax.plot(
        profile["r_mean_mm"],
        profile["pic_abs_mean_V_per_m"],
        color="#1f5aa6",
        lw=2.0,
        marker="o",
        ms=3.5,
        label="OPALX binned mean",
    )
    ax.plot(
        profile["r_mean_mm"],
        profile["ana_abs_mean_V_per_m"],
        color="#b24a3b",
        lw=2.0,
        marker="s",
        ms=3.2,
        label="analytic binned mean",
    )
    ax.set_xlabel("radius from diagnostic centroid [mm]")
    ax.set_ylabel(r"$|E|$ [V/m]")
    ax.set_title("Radial field magnitude")
    ax.grid(True, which="both", color="0.88", lw=0.7)
    ax.legend(frameon=True)
    fig.savefig(case_dir / "radial_field_vs_analytic.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ax.axhline(0.0, color="0.25", lw=1.0)
    ax.scatter(
        radius_mm,
        rel_mag[plot_index],
        s=7,
        alpha=0.14,
        color="#386cb0",
        linewidths=0,
        label="particles",
    )
    ax.plot(
        profile["r_mean_mm"],
        profile["rel_mag_mean"],
        color="#1f5aa6",
        lw=2.0,
        marker="o",
        ms=3.5,
        label="mean magnitude error",
    )
    ax.plot(
        profile["r_mean_mm"],
        profile["rel_mag_median"],
        color="#7b3294",
        lw=1.8,
        marker="s",
        ms=3.0,
        label="median magnitude error",
    )
    ax.set_xlabel("radius from diagnostic centroid [mm]")
    ax.set_ylabel(r"$(|E|_\mathrm{OPALX}-|E|_\mathrm{analytic})/|E|_\mathrm{analytic}$")
    ax.set_title("Radial field relative magnitude error")
    ax.grid(True, which="both", color="0.88", lw=0.7)
    ax.legend(frameon=True)
    fig.savefig(case_dir / "radial_field_relative_error.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ax.scatter(
        radius_mm,
        pic_er[plot_index],
        s=7,
        alpha=0.16,
        color="#386cb0",
        linewidths=0,
        label="OPALX particles",
    )
    ax.plot(
        profile["r_mean_mm"],
        profile["pic_er_mean_V_per_m"],
        color="#1f5aa6",
        lw=2.0,
        marker="o",
        ms=3.5,
        label="OPALX binned mean",
    )
    ax.plot(
        profile["r_mean_mm"],
        profile["ana_er_mean_V_per_m"],
        color="#b24a3b",
        lw=2.0,
        marker="s",
        ms=3.2,
        label="analytic binned mean",
    )
    ax.set_xlabel("radius from diagnostic centroid [mm]")
    ax.set_ylabel(r"$E_r$ [V/m]")
    ax.set_title("Signed radial field projection")
    ax.grid(True, which="both", color="0.88", lw=0.7)
    ax.legend(frameon=True)
    fig.savefig(case_dir / "signed_radial_field_vs_analytic.png", dpi=220)
    plt.close(fig)


def compute_case_metrics(
    raw: pd.DataFrame,
    comparator,
    sample_particles: int,
    sample_csv: Path,
    case_dir: Path,
) -> pd.DataFrame:
    step0 = raw.loc[raw["step"].eq(0) & raw["container"].eq(0)].copy()
    if step0.empty:
        raise RuntimeError("No c0 step-0 rows found in OPALX field dump")
    step0 = step0.sort_values(["rank", "local_index"]).reset_index(drop=True)
    center = (
        float(step0["x_m"].mean()),
        float(step0["y_m"].mean()),
        float(step0["z_m"].mean()),
    )

    pic = pd.DataFrame(
        {
            "x_m": step0["x_m"],
            "y_m": step0["y_m"],
            "z_m": step0["z_m"],
            "pic_Ex": step0["Ex_V_per_m"],
            "pic_Ey": step0["Ey_V_per_m"],
            "pic_Ez": step0["Ez_V_per_m"],
        }
    )
    pic["pic_E_abs"] = comparator.vector_abs(pic, "pic_")
    ax, ay, az = comparator.gaussian_e_field(
        pic["x_m"].to_numpy(float),
        pic["y_m"].to_numpy(float),
        pic["z_m"].to_numpy(float),
        charge=TOTAL_CHARGE_C,
        sigma=SIGMA_M,
        center=center,
    )
    pic["ana_Ex"] = ax
    pic["ana_Ey"] = ay
    pic["ana_Ez"] = az
    pic["ana_E_abs"] = comparator.vector_abs(pic, "ana_")
    comparison = comparator.add_error_columns(pic, rel_floor=1.0e-30, center=center)
    metrics = comparator.compute_metrics(comparison, rel_floor=1.0e-30)
    metrics["analytic_center_x_m"] = center[0]
    metrics["analytic_center_y_m"] = center[1]
    metrics["analytic_center_z_m"] = center[2]
    metrics["rms_x_m"] = float(step0["x_m"].std(ddof=0))
    metrics["rms_y_m"] = float(step0["y_m"].std(ddof=0))
    metrics["rms_z_m"] = float(step0["z_m"].std(ddof=0))

    save_case_radial_plots(comparison, center, case_dir, sample_particles)

    if sample_particles > 0:
        sample = comparison.sample(
            n=min(sample_particles, len(comparison)), random_state=20260703
        ).sort_index()
        sample.to_csv(sample_csv, index=False)
    return metrics


def save_aggregate_plots(summary: pd.DataFrame, output_dir: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output_dir / ".cache"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    grouped = summary.groupby(["n_grid", "nppg"], as_index=False).agg(
        relative_l2_mean=("relative_vector_L2_vs_analytic", "mean"),
        relative_l2_std=("relative_vector_L2_vs_analytic", "std"),
        n_particles=("configured_particles", "first"),
    )

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    for nppg, group in grouped.groupby("nppg"):
        ax.errorbar(
            group["n_grid"],
            group["relative_l2_mean"],
            yerr=group["relative_l2_std"].fillna(0.0),
            marker="o",
            capsize=3,
            label=f"NPPG={nppg}",
        )
    ref = grouped.groupby("n_grid", as_index=False)["relative_l2_mean"].mean().sort_values("n_grid")
    if not ref.empty:
        n0 = float(ref["n_grid"].iloc[-1])
        e0 = float(ref["relative_l2_mean"].iloc[-1])
        n = ref["n_grid"].to_numpy(float)
        ax.plot(n, e0 * (n0 / n) ** 2, "k--", lw=1.1, label=r"$O(N_\mathrm{grid}^{-2})$")
    ax.set_xlabel("N_GRID")
    ax.set_ylabel("relative vector L2 vs analytic")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="both", color="0.88", lw=0.6)
    ax.legend()
    fig.savefig(output_dir / "relative_l2_vs_grid.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    for n_grid, group in grouped.groupby("n_grid"):
        ax.errorbar(
            group["nppg"],
            group["relative_l2_mean"],
            yerr=group["relative_l2_std"].fillna(0.0),
            marker="o",
            capsize=3,
            label=f"N_GRID={n_grid}",
        )
    ax.set_xlabel("NPPG")
    ax.set_ylabel("relative vector L2 vs analytic")
    ax.set_yscale("log")
    ax.grid(True, which="both", color="0.88", lw=0.6)
    ax.legend()
    fig.savefig(output_dir / "relative_l2_vs_nppg.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ax.scatter(grouped["n_particles"], grouped["relative_l2_mean"], c=grouped["n_grid"], s=70)
    ax.set_xlabel("configured particles")
    ax.set_ylabel("mean relative vector L2 vs analytic")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", color="0.88", lw=0.6)
    fig.savefig(output_dir / "relative_l2_vs_particles.png", dpi=220)
    plt.close(fig)


def save_all_case_radial_plots(summary: pd.DataFrame, output_dir: Path) -> None:
    profiles = []
    for row in summary.itertuples(index=False):
        profile_path = output_dir / row.case / "radial_profile.csv"
        if not profile_path.exists():
            continue
        profile = pd.read_csv(profile_path)
        profile["case"] = row.case
        profile["n_grid"] = row.n_grid
        profile["nppg"] = row.nppg
        profile["seed"] = row.seed
        profiles.append(profile)
    if not profiles:
        return

    radial = pd.concat(profiles, ignore_index=True)
    n_grid_values = sorted(radial["n_grid"].unique())
    nppg_values = sorted(radial["nppg"].unique())

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    def make_grid_plot(
        filename: str,
        y_pic: str,
        y_ana: str | None,
        ylabel: str,
        title: str,
        zero_line: bool = False,
    ) -> None:
        fig, axes = plt.subplots(
            len(n_grid_values),
            len(nppg_values),
            figsize=(4.0 * len(nppg_values), 3.0 * len(n_grid_values)),
            sharex=True,
            constrained_layout=True,
        )
        axes = np.asarray(axes).reshape(len(n_grid_values), len(nppg_values))
        for i, n_grid in enumerate(n_grid_values):
            for j, nppg in enumerate(nppg_values):
                ax = axes[i, j]
                subset = radial[(radial["n_grid"] == n_grid) & (radial["nppg"] == nppg)]
                if subset.empty:
                    ax.set_axis_off()
                    continue
                if zero_line:
                    ax.axhline(0.0, color="0.30", lw=0.8)
                for seed, seed_group in subset.groupby("seed"):
                    seed_group = seed_group.sort_values("r_mean_mm")
                    ax.plot(
                        seed_group["r_mean_mm"],
                        seed_group[y_pic],
                        color="#1f5aa6",
                        lw=1.4,
                        alpha=0.35 if subset["seed"].nunique() > 1 else 1.0,
                        label="OPALX" if seed == subset["seed"].min() else None,
                    )
                if y_ana is not None:
                    first_seed = subset["seed"].min()
                    analytic = subset[subset["seed"] == first_seed].sort_values("r_mean_mm")
                    ax.plot(
                        analytic["r_mean_mm"],
                        analytic[y_ana],
                        color="#b24a3b",
                        lw=1.5,
                        ls="--",
                        label="analytic",
                    )
                ax.set_title(f"N_GRID={n_grid}, NPPG={nppg}")
                ax.grid(True, which="both", color="0.88", lw=0.6)
                if i == len(n_grid_values) - 1:
                    ax.set_xlabel("radius [mm]")
                if j == 0:
                    ax.set_ylabel(ylabel)
                if i == 0 and j == 0:
                    ax.legend(frameon=True, loc="best")
        fig.suptitle(title)
        fig.savefig(output_dir / filename, dpi=220)
        plt.close(fig)

    make_grid_plot(
        "radial_field_vs_analytic_all_cases.png",
        "pic_abs_mean_V_per_m",
        "ana_abs_mean_V_per_m",
        r"$|E|$ [V/m]",
        "Radial field magnitude by grid and particles per cell",
    )
    make_grid_plot(
        "signed_radial_field_vs_analytic_all_cases.png",
        "pic_er_mean_V_per_m",
        "ana_er_mean_V_per_m",
        r"$E_r$ [V/m]",
        "Signed radial field projection by grid and particles per cell",
    )
    make_grid_plot(
        "radial_field_relative_error_all_cases.png",
        "rel_mag_mean",
        None,
        r"mean relative $|E|$ error",
        "Radial field relative magnitude error by grid and particles per cell",
        zero_line=True,
    )


def main() -> None:
    args = parse_args()
    n_grid_values = parse_csv_ints(args.n_grid_values)
    nppg_values = parse_csv_ints(args.nppg_values)
    seeds = parse_csv_ints(args.seeds)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(args.output_dir / ".cache"))
    if args.output_dir.exists() and args.force:
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comparator = load_comparator_module()

    rows = []
    for n_grid in n_grid_values:
        for nppg in nppg_values:
            for seed in seeds:
                case = f"ngrid{n_grid:03d}_nppg{nppg:02d}_seed{seed}"
                case_dir = args.output_dir / case
                case_dir.mkdir(parents=True, exist_ok=True)
                deck = case_dir / "spacecharge_drift_30cm.in"
                configured_particles = render_input(n_grid, nppg, seed, deck)
                metrics_path = case_dir / "metrics.csv"
                if args.resume and metrics_path.exists():
                    print(f"Skipping {case}: metrics already exist", flush=True)
                    rows.append(pd.read_csv(metrics_path))
                    continue
                log = case_dir / "opalx.out"
                print(f"Running {case}: particles={configured_particles}", flush=True)
                run_opalx(args.opalx, deck, args.threads, log)
                raw = read_particle_field_h5(case_dir)
                metrics = compute_case_metrics(
                    raw,
                    comparator,
                    sample_particles=args.sample_particles,
                    sample_csv=case_dir / "particle_field_comparison_sample.csv",
                    case_dir=case_dir,
                )
                metrics.insert(0, "case", case)
                metrics.insert(1, "n_grid", n_grid)
                metrics.insert(2, "nppg", nppg)
                metrics.insert(3, "seed", seed)
                metrics.insert(4, "configured_particles", configured_particles)
                metrics["nominal_1_over_sqrt_N"] = 1.0 / math.sqrt(configured_particles)
                metrics["ratio_to_3_over_sqrt_N"] = (
                    metrics["relative_vector_L2_vs_analytic"]
                    / (3.0 / math.sqrt(configured_particles))
                )
                metrics.to_csv(metrics_path, index=False)
                rows.append(metrics)
                summary = pd.concat(rows, ignore_index=True)
                summary.to_csv(args.output_dir / "convergence_summary.csv", index=False)
                if not args.keep_raw:
                    for h5_file in case_dir.glob("*.h5"):
                        h5_file.unlink()

    summary = pd.concat(rows, ignore_index=True)
    summary.to_csv(args.output_dir / "convergence_summary.csv", index=False)
    grouped = summary.groupby(["n_grid", "nppg"], as_index=False).agg(
        configured_particles=("configured_particles", "first"),
        runs=("case", "count"),
        relative_l2_mean=("relative_vector_L2_vs_analytic", "mean"),
        relative_l2_std=("relative_vector_L2_vs_analytic", "std"),
        relative_l2_min=("relative_vector_L2_vs_analytic", "min"),
        relative_l2_max=("relative_vector_L2_vs_analytic", "max"),
        ratio_to_3_over_sqrt_N_mean=("ratio_to_3_over_sqrt_N", "mean"),
    )
    grouped.to_csv(args.output_dir / "convergence_grouped_summary.csv", index=False)
    save_aggregate_plots(summary, args.output_dir)
    save_all_case_radial_plots(summary, args.output_dir)
    (args.output_dir / "README.md").write_text(
        "\n".join(
            [
                "# Drift Space-Charge Convergence Study",
                "",
                f"N_GRID values: `{list(n_grid_values)}`",
                f"NPPG values: `{list(nppg_values)}`",
                f"Seeds: `{list(seeds)}`",
                "",
                "Files:",
                "",
                "- `convergence_summary.csv`: one row per OPALX run.",
                "- `convergence_grouped_summary.csv`: grouped mean/std by `(N_GRID, NPPG)`.",
                "- `relative_l2_vs_grid.png`, `relative_l2_vs_nppg.png`,",
                "  `relative_l2_vs_particles.png`: aggregate convergence plots.",
                "- `radial_field_vs_analytic_all_cases.png`,",
                "  `signed_radial_field_vs_analytic_all_cases.png`, and",
                "  `radial_field_relative_error_all_cases.png`: radial profile grids.",
                "- Per-case directories contain `opalx.out`, the rendered input deck,",
                "  `metrics.csv`, `radial_profile.csv`, radial PNGs, and a deterministic",
                "  sampled comparison CSV.",
                "",
                "Raw H5 particle dumps are deleted after each successful case unless",
                "`--keep-raw` is used.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
