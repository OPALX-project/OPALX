#!/usr/bin/env python3
"""Analyze the staged BeamBeam witness timing smoke run.

The script is intentionally lightweight: it reads the H5 particle/field dumps
and the ASCII SDDS stat files produced by ``spacecharge_drift_withness.in`` and
creates timing/mesh plots that can be regenerated after deck or branch changes.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CLIGHT = 299_792_458.0
EMASS_GEV = 0.00051099895


def read_h5_summary(path: Path) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    with h5py.File(path, "r") as h5:
        for key in sorted(h5.keys(), key=lambda name: int(name.split("#")[1])):
            group = h5[key]
            attrs = group.attrs
            n = len(group["x"])
            row = {
                "h5_step": int(attrs["Step"][0]),
                "global_step": int(attrs["GlobalTrackStep"][0]),
                "time_s": float(attrs["TIME"][0]),
                "time_ps": float(attrs["TIME"][0]) * 1.0e12,
                "spos_m": float(attrs["SPOS"][0]),
                "num_particles": n,
                "charge_c": float(attrs["CHARGE"][0]),
                "rms_x_m": float(attrs["RMSX"][0]),
                "rms_y_m": float(attrs["RMSX"][1]),
                "rms_s_m": float(attrs["RMSX"][2]),
                "min_x_m": float(attrs["minX"][0]),
                "min_y_m": float(attrs["minX"][1]),
                "min_s_m": float(attrs["minX"][2]),
                "max_x_m": float(attrs["maxX"][0]),
                "max_y_m": float(attrs["maxX"][1]),
                "max_s_m": float(attrs["maxX"][2]),
                "ref_z_m": float(attrs["RefPartR"][2]),
                "mean_x_m": float(attrs["centroid"][0]),
                "mean_y_m": float(attrs["centroid"][1]),
                "mean_s_m": float(attrs["centroid"][2]),
            }
            for name in ("x", "y", "z", "px", "py", "pz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"):
                if name in group and n:
                    values = group[name][...]
                    row[f"{name}_0"] = float(values[0])
                    row[f"{name}_mean"] = float(np.mean(values))
            rows.append(row)
    return pd.DataFrame(rows)


def read_stat(path: Path) -> pd.DataFrame:
    text = path.read_text()
    columns = re.findall(r"&column\s+name=([^,\n]+),", text)
    data_start = text.index("&data")
    lines = [line.strip() for line in text[data_start:].splitlines() if line.strip()]
    end_idx = lines.index("&end")
    payload = lines[end_idx + 1 :]
    if len(payload) < 3:
        return pd.DataFrame(columns=columns)
    rows = []
    for line in payload[3:]:
        values = line.split()
        if len(values) != len(columns):
            continue
        rows.append(values)
    df = pd.DataFrame(rows, columns=columns)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def primary_beta(primary_kinetic_energy_gev: float) -> float:
    gamma = (primary_kinetic_energy_gev + EMASS_GEV) / EMASS_GEV
    return math.sqrt(1.0 - 1.0 / (gamma * gamma))


def infer_primary_source_r0z(
    run_dir: Path,
    beta_primary: float,
    bb_ip_s: float,
    witness_t0: float,
    witness_ct_m: float,
    primary_sigma_z_m: float,
    explicit_r0z: float | None,
) -> float:
    if explicit_r0z is not None:
        return explicit_r0z

    deck = run_dir / "spacecharge_drift_withness.in"
    if deck.exists():
        text = deck.read_text()
        if (
            "primary_centroid_at_witness_t0" in text
            and "track12_edge_reference_s" in text
            and "primary_edge_sigma" in text
        ):
            edge_reference_s = bb_ip_s - 1.5 * primary_sigma_z_m
            centroid_at_t0 = edge_reference_s - 3.0 * primary_sigma_z_m
            return centroid_at_t0 - beta_primary * CLIGHT * witness_t0
        if "primary_source_r0z" in text and "R0Z = primary_source_r0z" in text:
            return bb_ip_s - beta_primary * CLIGHT * witness_t0 - witness_ct_m

    return 0.0


def nearest_row(df: pd.DataFrame, column: str, value: float) -> pd.Series:
    idx = (df[column] - value).abs().idxmin()
    return df.loc[idx]


def plot_timing(
    c0: pd.DataFrame,
    c1_stat: pd.DataFrame,
    c2_stat: pd.DataFrame,
    c1_h5: pd.DataFrame,
    out: Path,
    bb_ip_s: float,
    witness_t0: float,
    retire_time: float | None,
    primary_sigma_z_m: float,
    track12_edge_reference_s: float,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.4), sharex=True)
    ax = axes[0]
    c0_s_column = "absolute_spos_m" if "absolute_spos_m" in c0.columns else "spos_m"
    c0_centroid_s = c0[c0_s_column]
    ax.plot(c0["time_ps"], 1.0e3 * c0_centroid_s, "o-", ms=4, lw=1.5, label="c0 centroid H5")
    ax.plot(
        c0["time_ps"],
        1.0e3 * (c0_centroid_s + 3.0 * primary_sigma_z_m),
        "o--",
        ms=3,
        lw=1.1,
        label="c0 +3sigma_z edge",
    )
    ax.axhline(1.0e3 * bb_ip_s, color="black", ls="--", lw=1.1, label="BeamBeam IP")
    ax.axhline(
        1.0e3 * track12_edge_reference_s,
        color="tab:orange",
        ls=":",
        lw=1.4,
        label="track12 first witness",
    )
    ax.axvline(1.0e12 * witness_t0, color="tab:red", ls="--", lw=1.1, label="witness T0")
    if retire_time is not None:
        ax.axvline(1.0e12 * retire_time, color="tab:green", ls="--", lw=1.1, label="c0 retire")
    first_witness = c1_h5.iloc[0]
    ax.axvline(first_witness["time_ps"], color="tab:purple", ls=":", lw=1.3, label="first c1/c2 dump")
    ax.set_ylabel("c0 $s$ [mm]")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ax.step(c1_stat["t"] * 1.0e3, c1_stat["numParticles"], where="post", lw=1.6, label="c1 stat")
    ax.step(c2_stat["t"] * 1.0e3, c2_stat["numParticles"], where="post", lw=1.2, ls="--", label="c2 stat")
    ax.axvline(1.0e12 * witness_t0, color="tab:red", ls="--", lw=1.1)
    if retire_time is not None:
        ax.axvline(1.0e12 * retire_time, color="tab:green", ls="--", lw=1.1)
    ax.set_ylabel("witness macroparticles")
    max_particles = max(float(c1_stat["numParticles"].max()), float(c2_stat["numParticles"].max()), 1.0)
    ax.set_ylim(-0.05 * max_particles, 1.15 * max_particles)
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    if {"Ex_0", "Ey_0", "Ez_0"}.issubset(c1_h5.columns):
        e_abs = np.sqrt(c1_h5["Ex_0"] ** 2 + c1_h5["Ey_0"] ** 2 + c1_h5["Ez_0"] ** 2)
        ax.plot(c1_h5["time_ps"], e_abs, "o-", ms=4, lw=1.5, label="c1 |E| at particle")
        ax.plot(c1_h5["time_ps"], c1_h5["Ex_0"], "s-", ms=3, lw=1.0, label="c1 Ex")
        ax.plot(c1_h5["time_ps"], c1_h5["Ey_0"], "^-", ms=3, lw=1.0, label="c1 Ey")
    ax.axvline(1.0e12 * witness_t0, color="tab:red", ls="--", lw=1.1)
    if retire_time is not None:
        ax.axvline(1.0e12 * retire_time, color="tab:green", ls="--", lw=1.1)
    ax.set_ylabel("field [V/m]")
    ax.set_xlabel("time [ps]")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.legend(loc="best", fontsize=8, ncol=3)
    ax.grid(True, alpha=0.25)

    fig.suptitle("spacecharge_drift_withness timing")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_mesh(
    c0_h5_path: Path,
    c0_summary: pd.DataFrame,
    out: Path,
    c0_mesh_row: pd.Series,
    aperture_radius: float,
    nxy: int,
) -> None:
    row = c0_mesh_row
    c0_s_m = row["absolute_spos_m"] if "absolute_spos_m" in c0_summary.columns else row["spos_m"]
    h5_key = f"Step#{int(row['h5_step'])}"
    with h5py.File(c0_h5_path, "r") as h5:
        group = h5[h5_key]
        x = group["x"][...]
        y = group["y"][...]

    dx = 2.0 * aperture_radius / nxy
    edges = np.linspace(-aperture_radius, aperture_radius, nxy + 1)
    theta = np.linspace(0.0, 2.0 * np.pi, 512)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    ax = axes[0]
    ax.plot(aperture_radius * np.cos(theta) * 1.0e6, aperture_radius * np.sin(theta) * 1.0e6, color="black", lw=1.0)
    for edge in edges:
        ax.axvline(edge * 1.0e6, color="0.85", lw=0.4, zorder=0)
        ax.axhline(edge * 1.0e6, color="0.85", lw=0.4, zorder=0)
    ax.scatter(x * 1.0e6, y * 1.0e6, s=16, color="tab:blue", alpha=0.75)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_title(f"full aperture, dx={dx * 1.0e6:.2f} um")
    ax.grid(False)

    ax = axes[1]
    zoom = 15.0e-6
    zoom_edges = edges[(edges >= -zoom) & (edges <= zoom)]
    for edge in zoom_edges:
        ax.axvline(edge * 1.0e6, color="0.80", lw=0.8, zorder=0)
        ax.axhline(edge * 1.0e6, color="0.80", lw=0.8, zorder=0)
    ax.scatter(x * 1.0e6, y * 1.0e6, s=24, color="tab:blue", alpha=0.78, label="c0 macroparticles")
    rms_x = row["rms_x_m"] * 1.0e6
    rms_y = row["rms_y_m"] * 1.0e6
    ellipse = patches.Ellipse((0.0, 0.0), 2.0 * rms_x, 2.0 * rms_y, fill=False, color="tab:red", lw=1.4, label="1 rms")
    ax.add_patch(ellipse)
    ax.set_xlim(-zoom * 1.0e6, zoom * 1.0e6)
    ax.set_ylim(-zoom * 1.0e6, zoom * 1.0e6)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_title("injection zoom")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(False)

    fig.suptitle(
        f"nearest saved c0 transverse mesh at t={row['time_ps']:.3f} ps, "
        f"s={c0_s_m * 1.0e3:.4f} mm"
    )
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def momentum_to_beta(px: float, py: float, pz: float) -> tuple[float, float, float, float]:
    p2 = px * px + py * py + pz * pz
    gamma = math.sqrt(1.0 + p2)
    return px / gamma, py / gamma, pz / gamma, gamma


def make_witness_kinematics_summary(c1: pd.DataFrame, c2: pd.DataFrame, aperture_radius: float) -> pd.DataFrame:
    rows = []
    for name, df in (("c1", c1), ("c2", c2)):
        first = df.iloc[0]
        last = df.iloc[-1]
        beta_x, beta_y, beta_z, gamma = momentum_to_beta(first["px_0"], first["py_0"], first["pz_0"])
        rows.append(
            {
                "container": name,
                "first_time_ps": first["time_ps"],
                "first_global_step": int(first["global_step"]),
                "first_x_um": first["x_0"] * 1.0e6,
                "first_y_um": first["y_0"] * 1.0e6,
                "first_z_um": first["z_0"] * 1.0e6,
                "first_px_norm": first["px_0"],
                "first_py_norm": first["py_0"],
                "first_pz_norm": first["pz_0"],
                "first_gamma_from_p": gamma,
                "first_beta_x": beta_x,
                "first_beta_y": beta_y,
                "first_beta_z": beta_z,
                "final_time_ps": last["time_ps"],
                "final_global_step": int(last["global_step"]),
                "final_x_um": last["x_0"] * 1.0e6,
                "final_y_um": last["y_0"] * 1.0e6,
                "final_z_um": last["z_0"] * 1.0e6,
                "max_abs_x_um": max(abs(df["x_0"].min()), abs(df["x_0"].max())) * 1.0e6,
                "max_abs_y_um": max(abs(df["y_0"].min()), abs(df["y_0"].max())) * 1.0e6,
                "max_abs_y_over_aperture": max(abs(df["y_0"].min()), abs(df["y_0"].max())) / aperture_radius,
            }
        )
    return pd.DataFrame(rows)


def plot_witness_motion(c1: pd.DataFrame, c2: pd.DataFrame, out: Path, aperture_radius: float) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))

    ax = axes[0]
    theta = np.linspace(0.0, 2.0 * np.pi, 512)
    ax.plot(
        aperture_radius * np.cos(theta) * 1.0e6,
        aperture_radius * np.sin(theta) * 1.0e6,
        color="black",
        ls="--",
        lw=1.0,
        label="c0 source aperture",
    )
    ax.plot(c1["x_0"] * 1.0e6, c1["y_0"] * 1.0e6, "o-", ms=4, lw=1.4, label="c1/e-")
    ax.plot(c2["x_0"] * 1.0e6, c2["y_0"] * 1.0e6, "s-", ms=4, lw=1.4, label="c2/e+")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_title("transverse path")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    t0 = min(c1["time_ps"].iloc[0], c2["time_ps"].iloc[0])
    ax.plot(c1["time_ps"] - t0, c1["y_0"] * 1.0e6, "o-", ms=4, lw=1.4, label="c1/e- y")
    ax.plot(c2["time_ps"] - t0, c2["y_0"] * 1.0e6, "s-", ms=4, lw=1.4, label="c2/e+ y")
    ax.axhline(aperture_radius * 1.0e6, color="black", ls="--", lw=1.0)
    ax.axhline(-aperture_radius * 1.0e6, color="black", ls="--", lw=1.0)
    ax.set_xlabel("time since first witness dump [ps]")
    ax.set_ylabel("y [um]")
    ax.set_title("motion relative to compact source aperture")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.25)

    fig.suptitle("transverse witness motion")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def line_gaussian_relative_field(
    time_since_t0_s: np.ndarray, beta_witness: float, beta_primary: float, b_m: float, sigma_z_m: float
) -> np.ndarray:
    t = np.maximum(time_since_t0_s, 0.0)
    rho = np.sqrt(b_m * b_m + (beta_witness * CLIGHT * t) ** 2)
    z = beta_primary * CLIGHT * t
    return (b_m / rho) * np.exp(-0.5 * (z / sigma_z_m) ** 2)


def plot_near_ip_cutoff(
    c1: pd.DataFrame,
    c2: pd.DataFrame,
    out_png: Path,
    out_csv: Path,
    witness_t0: float,
    retire_time: float,
    beta_witness: float,
    beta_primary: float,
    b_m: float,
    sigma_z_m: float,
) -> None:
    rows = []
    for name, df in (("c1", c1), ("c2", c2)):
        e_abs = np.sqrt(df["Ex_0"] ** 2 + df["Ey_0"] ** 2 + df["Ez_0"] ** 2)
        e0 = e_abs.iloc[0] if e_abs.iloc[0] != 0.0 else np.nan
        for _, row in df.iterrows():
            t_since = row["time_s"] - witness_t0
            e_value = math.sqrt(row["Ex_0"] ** 2 + row["Ey_0"] ** 2 + row["Ez_0"] ** 2)
            estimate = line_gaussian_relative_field(
                np.asarray([t_since]), beta_witness, beta_primary, b_m, sigma_z_m
            )[0]
            rows.append(
                {
                    "container": name,
                    "time_ps": row["time_ps"],
                    "time_since_t0_fs": t_since * 1.0e15,
                    "E_abs_V_per_m": e_value,
                    "E_abs_relative_to_first": e_value / e0,
                    "line_gaussian_relative_estimate": estimate,
                    "retire_time_ps": retire_time * 1.0e12,
                }
            )
    summary = pd.DataFrame(rows)
    summary.to_csv(out_csv, index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for name, marker in (("c1", "o"), ("c2", "s")):
        subset = summary[summary["container"] == name]
        ax.semilogy(
            subset["time_since_t0_fs"],
            subset["E_abs_relative_to_first"],
            marker + "-",
            ms=4,
            lw=1.2,
            label=f"{name} OPALX |E| / first",
        )

    max_time_s = max(1.0e-15, summary["time_since_t0_fs"].max() * 1.0e-15)
    t_grid = np.linspace(0.0, max_time_s, 512)
    ax.semilogy(
        t_grid * 1.0e15,
        line_gaussian_relative_field(t_grid, beta_witness, beta_primary, b_m, sigma_z_m),
        color="black",
        ls="--",
        lw=1.4,
        label="line-Gaussian estimate",
    )
    ax.axhline(0.01, color="0.35", ls=":", lw=1.2, label="1%")
    ax.axvline((retire_time - witness_t0) * 1.0e15, color="tab:green", ls="--", lw=1.2, label="c0 retire")
    ax.set_xlabel("time since witness T0 [fs]")
    ax.set_ylabel("relative field")
    ax.set_title("near-IP active-field cutoff")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path("sandbox/Drift-Experiment/withness_timing_pair4_staged_fixed"),
    )
    parser.add_argument("--bb-ip-s", type=float, default=4.0e-3)
    parser.add_argument("--primary-kinetic-energy-gev", type=float, default=0.245)
    parser.add_argument("--witness-kinetic-energy-gev", type=float, default=EMASS_GEV)
    parser.add_argument("--primary-sigma-xy-m", type=float, default=1.944325075701e-6)
    parser.add_argument("--primary-sigma-z-m", type=float, default=0.6e-3)
    parser.add_argument("--witness-ct-m", type=float, default=0.0)
    parser.add_argument("--witness-t0-s", type=float, default=None)
    parser.add_argument("--primary-source-r0z-m", type=float, default=None)
    parser.add_argument("--near-ip-active-time-s", type=float, default=1.6e-12)
    parser.add_argument("--aperture-radius-m", type=float, default=5.0e-4)
    parser.add_argument("--nxy", type=int, default=16)
    args = parser.parse_args()

    run_dir = args.run_dir
    beta = primary_beta(args.primary_kinetic_energy_gev)
    witness_t0 = (
        args.witness_t0_s
        if args.witness_t0_s is not None
        else (args.bb_ip_s - args.witness_ct_m) / (beta * CLIGHT)
    )
    retire_time = witness_t0 + args.near_ip_active_time_s
    witness_beta = primary_beta(args.witness_kinetic_energy_gev)
    primary_source_r0z = infer_primary_source_r0z(
        run_dir,
        beta,
        args.bb_ip_s,
        witness_t0,
        args.witness_ct_m,
        args.primary_sigma_z_m,
        args.primary_source_r0z_m,
    )
    track12_edge_reference_s = args.bb_ip_s - 1.5 * args.primary_sigma_z_m

    c0 = read_h5_summary(run_dir / "spacecharge_drift_withness_c0.h5")
    c0["absolute_spos_m"] = c0["spos_m"] + primary_source_r0z
    c1 = read_h5_summary(run_dir / "spacecharge_drift_withness_c1.h5")
    c2 = read_h5_summary(run_dir / "spacecharge_drift_withness_c2.h5")
    c1_stat = read_stat(run_dir / "spacecharge_drift_withness_c1.stat")
    c2_stat = read_stat(run_dir / "spacecharge_drift_withness_c2.stat")

    first_witness = c1.iloc[0]
    c0_at_first_witness = nearest_row(c0, "global_step", first_witness["global_step"])
    c0_at_ip = nearest_row(c0, "absolute_spos_m", args.bb_ip_s)

    timing_summary = pd.DataFrame(
        [
            {
                "witness_t0_ps": witness_t0 * 1.0e12,
                "first_witness_h5_time_ps": first_witness["time_ps"],
                "first_witness_global_step": int(first_witness["global_step"]),
                "first_witness_minus_t0_fs": (first_witness["time_s"] - witness_t0) * 1.0e15,
                "nearest_c0_to_first_witness_global_step": int(c0_at_first_witness["global_step"]),
                "nearest_c0_to_first_witness_delta_steps": int(c0_at_first_witness["global_step"] - first_witness["global_step"]),
                "nearest_c0_to_first_witness_delta_fs": (c0_at_first_witness["time_s"] - first_witness["time_s"]) * 1.0e15,
                "configured_retire_time_ps": retire_time * 1.0e12,
                "configured_retire_minus_t0_fs": args.near_ip_active_time_s * 1.0e15,
                "primary_source_r0z_mm": primary_source_r0z * 1.0e3,
                "track12_edge_reference_s_mm": track12_edge_reference_s * 1.0e3,
                "c0_centroid_at_witness_t0_mm": (primary_source_r0z + beta * CLIGHT * witness_t0) * 1.0e3,
                "c0_plus_3sigma_edge_at_witness_t0_mm": (
                    primary_source_r0z + beta * CLIGHT * witness_t0 + 3.0 * args.primary_sigma_z_m
                )
                * 1.0e3,
                "c0_s_at_first_witness_mm": c0_at_first_witness["absolute_spos_m"] * 1.0e3,
                "c0_s_minus_ip_at_first_witness_um": (c0_at_first_witness["absolute_spos_m"] - args.bb_ip_s) * 1.0e6,
                "nearest_c0_ip_time_ps": c0_at_ip["time_ps"],
                "nearest_c0_ip_s_mm": c0_at_ip["absolute_spos_m"] * 1.0e3,
                "nearest_c0_ip_global_step": int(c0_at_ip["global_step"]),
                "c0_rms_x_um_at_injection": c0_at_first_witness["rms_x_m"] * 1.0e6,
                "c0_rms_y_um_at_injection": c0_at_first_witness["rms_y_m"] * 1.0e6,
                "c0_rms_s_mm_at_injection": c0_at_first_witness["rms_s_m"] * 1.0e3,
                "c0_max_abs_x_um_at_injection": max(abs(c0_at_first_witness["min_x_m"]), abs(c0_at_first_witness["max_x_m"])) * 1.0e6,
                "c0_max_abs_y_um_at_injection": max(abs(c0_at_first_witness["min_y_m"]), abs(c0_at_first_witness["max_y_m"])) * 1.0e6,
                "c0_max_abs_s_mm_at_injection": max(abs(c0_at_first_witness["min_s_m"]), abs(c0_at_first_witness["max_s_m"])) * 1.0e3,
                "aperture_radius_um": args.aperture_radius_m * 1.0e6,
                "nxy": args.nxy,
                "transverse_cell_size_um": 2.0 * args.aperture_radius_m / args.nxy * 1.0e6,
            }
        ]
    )
    timing_summary.to_csv(run_dir / "timing_mesh_summary.csv", index=False)

    field_rows = c1.merge(
        c2[["global_step", "Ex_0", "Ey_0", "Ez_0", "Bx_0", "By_0", "Bz_0"]],
        on="global_step",
        suffixes=("_c1", "_c2"),
    )
    field_rows["E_abs_c1"] = np.sqrt(field_rows["Ex_0_c1"] ** 2 + field_rows["Ey_0_c1"] ** 2 + field_rows["Ez_0_c1"] ** 2)
    field_rows["E_abs_c2"] = np.sqrt(field_rows["Ex_0_c2"] ** 2 + field_rows["Ey_0_c2"] ** 2 + field_rows["Ez_0_c2"] ** 2)
    field_rows.to_csv(run_dir / "witness_field_samples.csv", index=False)

    witness_kinematics = make_witness_kinematics_summary(c1, c2, args.aperture_radius_m)
    witness_kinematics.to_csv(run_dir / "witness_kinematics_summary.csv", index=False)

    plot_timing(
        c0,
        c1_stat,
        c2_stat,
        c1,
        run_dir / "witness_timing_overview.png",
        args.bb_ip_s,
        witness_t0,
        retire_time,
        args.primary_sigma_z_m,
        track12_edge_reference_s,
    )
    plot_mesh(
        run_dir / "spacecharge_drift_withness_c0.h5",
        c0,
        run_dir / "c0_injection_mesh_xy.png",
        c0_at_first_witness,
        args.aperture_radius_m,
        args.nxy,
    )
    plot_witness_motion(c1, c2, run_dir / "witness_transverse_motion.png", args.aperture_radius_m)
    plot_near_ip_cutoff(
        c1,
        c2,
        run_dir / "near_ip_field_cutoff.png",
        run_dir / "near_ip_field_cutoff.csv",
        witness_t0,
        retire_time,
        witness_beta,
        beta,
        args.primary_sigma_xy_m,
        args.primary_sigma_z_m,
    )

    print(timing_summary.to_string(index=False))
    print(witness_kinematics.to_string(index=False))
    print(f"Wrote {run_dir / 'witness_timing_overview.png'}")
    print(f"Wrote {run_dir / 'c0_injection_mesh_xy.png'}")
    print(f"Wrote {run_dir / 'witness_transverse_motion.png'}")
    print(f"Wrote {run_dir / 'near_ip_field_cutoff.png'}")


if __name__ == "__main__":
    main()
