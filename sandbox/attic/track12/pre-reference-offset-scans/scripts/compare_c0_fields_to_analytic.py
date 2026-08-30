#!/usr/bin/env python3
"""Compare OPALX c0 primary fields against the smooth manufactured Gaussian."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TRACK12_MODULE = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_OUTPUT_DIR = ROOT / "sandbox" / "track12particles" / "opalx" / "c0_field_compare"


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TRACK12_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure_matplotlib(output_dir: Path) -> None:
    cache_root = output_dir / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def source_for_pair4(module):
    sources = module.make_anisotropic_sources(
        module.DEFAULT_SOURCE_KINETIC_MEV,
        module.DEFAULT_SOURCE_CHARGE_C,
        module.SIGMA_XY_M,
        module.SIGMA_XY_M,
        module.SIGMA_Z_M,
        "copropagating",
        0.0,
        module.DEFAULT_MC_SOURCE_PARTICLES,
        module.DEFAULT_MC_SEED,
    )
    # Pair 4 has ct=0, so no additional source shift is needed.
    return sources[0]


def analytic_fields(rows: pd.DataFrame, module, source) -> pd.DataFrame:
    values = []
    for row in rows.itertuples(index=False):
        position = np.array([row.x_m, row.y_m, row.s_minus_ip_m], dtype=float)
        e_field, b_field = module.anisotropic_source_lab_fields(position, 0.0, source)
        values.append(
            {
                "Ex_analytic": e_field[0],
                "Ey_analytic": e_field[1],
                "Ez_analytic": e_field[2],
                "Bx_analytic": b_field[0],
                "By_analytic": b_field[1],
                "Bz_analytic": b_field[2],
                "E_abs_analytic": float(np.linalg.norm(e_field)),
                "B_abs_analytic": float(np.linalg.norm(b_field)),
            }
        )
    return pd.DataFrame(values)


def summarize(data: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for column in [
        "E_abs_opalx",
        "E_abs_analytic",
        "E_abs_ratio",
        "B_abs_opalx",
        "B_abs_analytic",
        "B_abs_ratio",
        "E_direction_cosine",
        "B_direction_cosine",
        "Ex_opalx",
        "Ex_analytic",
        "By_opalx",
        "By_analytic",
    ]:
        values = data[column].replace([np.inf, -np.inf], np.nan).dropna()
        rows.append(
            {
                "case": label,
                "quantity": column,
                "count": int(values.size),
                "median": values.median(),
                "p05": values.quantile(0.05),
                "p95": values.quantile(0.95),
                "mean": values.mean(),
            }
        )
    return pd.DataFrame(rows)


def analyze_case(path: Path, label: str, sample: int, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    module = load_track12_module()
    source = source_for_pair4(module)
    data = pd.read_csv(path)
    data = data.loc[data["container"].eq(0) & data["step"].eq(0)].copy()
    if sample > 0 and len(data) > sample:
        data = data.sample(n=sample, random_state=20260702).sort_values("id").reset_index(drop=True)
    else:
        data = data.reset_index(drop=True)

    analytic = analytic_fields(data, module, source)
    out = pd.concat([data, analytic], axis=1)
    out["E_abs_ratio"] = out["E_abs_V_per_m"] / out["E_abs_analytic"]
    out["B_abs_ratio"] = out["B_abs_T"] / out["B_abs_analytic"]
    e_dot = (
        out["Ex_V_per_m"] * out["Ex_analytic"]
        + out["Ey_V_per_m"] * out["Ey_analytic"]
        + out["Ez_V_per_m"] * out["Ez_analytic"]
    )
    b_dot = (
        out["Bx_T"] * out["Bx_analytic"]
        + out["By_T"] * out["By_analytic"]
        + out["Bz_T"] * out["Bz_analytic"]
    )
    out["E_direction_cosine"] = e_dot / (out["E_abs_V_per_m"] * out["E_abs_analytic"])
    out["B_direction_cosine"] = b_dot / (out["B_abs_T"] * out["B_abs_analytic"])
    out = out.rename(
        columns={
            "Ex_V_per_m": "Ex_opalx",
            "Ey_V_per_m": "Ey_opalx",
            "Ez_V_per_m": "Ez_opalx",
            "Bx_T": "Bx_opalx",
            "By_T": "By_opalx",
            "Bz_T": "Bz_opalx",
            "E_abs_V_per_m": "E_abs_opalx",
            "B_abs_T": "B_abs_opalx",
        }
    )
    out["case"] = label
    out_path = output_dir / f"{label}_c0_field_sample.csv"
    out.to_csv(out_path, index=False)
    return out, summarize(out, label)


def plot(combined: pd.DataFrame, output_dir: Path) -> None:
    configure_matplotlib(output_dir)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.6, 4.6), constrained_layout=True)
    for label, group in combined.groupby("case", sort=False):
        ax.hist(
            np.log10(group["E_abs_ratio"].replace([np.inf, -np.inf], np.nan).dropna()),
            bins=80,
            histtype="step",
            lw=1.2,
            label=label,
        )
    ax.axvline(0.0, color="0.2", ls="--", lw=1.0)
    ax.set_xlabel(r"$\log_{10}(|E|_\mathrm{OPALX}/|E|_\mathrm{analytic})$")
    ax.set_ylabel("sampled source particles")
    ax.grid(True, color="0.88", lw=0.6)
    ax.legend(fontsize=8)
    fig.savefig(output_dir / "c0_e_abs_ratio_hist.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 4.6), constrained_layout=True)
    for label, group in combined.groupby("case", sort=False):
        ax.hist(
            group["E_direction_cosine"].replace([np.inf, -np.inf], np.nan).dropna(),
            bins=80,
            histtype="step",
            lw=1.2,
            label=label,
        )
    ax.axvline(1.0, color="0.2", ls="--", lw=1.0)
    ax.set_xlabel("E direction cosine, OPALX vs analytic")
    ax.set_ylabel("sampled source particles")
    ax.grid(True, color="0.88", lw=0.6)
    ax.legend(fontsize=8)
    fig.savefig(output_dir / "c0_e_direction_cosine_hist.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-csv", type=Path, required=True)
    parser.add_argument("--free-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sample", type=int, default=20000)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fixed, fixed_summary = analyze_case(args.fixed_csv, "fixed_window", args.sample, args.output_dir)
    free, free_summary = analyze_case(args.free_csv, "default_following_mesh", args.sample, args.output_dir)
    combined = pd.concat([fixed, free], ignore_index=True)
    combined.to_csv(args.output_dir / "c0_field_sample_combined.csv", index=False)
    summary = pd.concat([fixed_summary, free_summary], ignore_index=True)
    summary.to_csv(args.output_dir / "c0_field_summary.csv", index=False)
    plot(combined, args.output_dir)
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.6g}"))
    print(f"wrote {args.output_dir}")


if __name__ == "__main__":
    main()
