#!/usr/bin/env python3
"""Analyze gamma-gamma pair FROMFILE momentum directions.

The FROMFILE columns are interpreted as

    x y z px py pz

where px, py, pz are beta*gamma = p/(m_e c).  Files may also contain extra
columns such as t.  When t is present, the script also prints time diagnostics
and can write a time histogram.  By default the script prints numeric direction
statistics only.  Use --plot to also write diagnostic figures; plotting is kept
optional because some macOS/font-cache environments can fail while initializing
matplotlib.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


ME_KEV = 510.99895
COLUMNS = ["x", "y", "z", "px", "py", "pz"]
TIME_COLUMN = "t"


def read_fromfile(path: Path) -> pd.DataFrame:
    """Read an OPALX FROMFILE distribution into a dataframe."""
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} is too short to be a FROMFILE distribution")

    n_declared = int(lines[0].strip())
    header = lines[1].split()
    missing = [column for column in COLUMNS if column not in header]
    if missing:
        raise ValueError(f"{path}: missing required FROMFILE columns {missing}")
    df = pd.read_csv(path, sep=r"\s+", skiprows=2, names=header)
    if len(df) != n_declared:
        raise ValueError(f"{path}: declared {n_declared} particles, read {len(df)}")
    output_columns = COLUMNS + ([TIME_COLUMN] if TIME_COLUMN in header else [])
    return df.loc[:, output_columns].copy()


def add_momentum_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add beta-gamma magnitude, unit direction, angles, and kinetic energy."""
    out = df.copy()
    p = out[["px", "py", "pz"]].to_numpy(float)
    p_abs = np.linalg.norm(p, axis=1)
    direction = p / p_abs[:, None]

    out["p_abs"] = p_abs
    out["pt"] = np.linalg.norm(p[:, :2], axis=1)
    out["ux"] = direction[:, 0]
    out["uy"] = direction[:, 1]
    out["uz"] = direction[:, 2]
    out["theta_deg"] = np.degrees(np.arccos(np.clip(out["uz"], -1.0, 1.0)))
    out["phi_deg"] = np.degrees(np.arctan2(out["uy"], out["ux"]))

    gamma = np.sqrt(1.0 + out["p_abs"] ** 2)
    out["kinetic_keV"] = (gamma - 1.0) * ME_KEV
    return out


def circular_mean_deg(phi_deg: pd.Series) -> float:
    """Return circular mean angle in degrees."""
    phi = np.radians(phi_deg.to_numpy(float))
    return float(np.degrees(np.arctan2(np.sin(phi).mean(), np.cos(phi).mean())))


def circular_resultant(phi_deg: pd.Series) -> float:
    """Return |<exp(i phi)>| as an azimuthal anisotropy measure."""
    phi = np.radians(phi_deg.to_numpy(float))
    return float(np.hypot(np.cos(phi).mean(), np.sin(phi).mean()))


def summarize_species(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build a compact per-species summary table."""
    rows = []
    for species, df in frames.items():
        u_mean = df[["ux", "uy", "uz"]].mean()
        p_mean = df[["px", "py", "pz"]].mean()
        values = {
            "N": len(df),
            "<px> beta_gamma": p_mean["px"],
            "<py> beta_gamma": p_mean["py"],
            "<pz> beta_gamma": p_mean["pz"],
            "<ux>": u_mean["ux"],
            "<uy>": u_mean["uy"],
            "<uz>": u_mean["uz"],
            "|<u>|": np.linalg.norm(u_mean.to_numpy()),
            "pz>0 frac": (df["pz"] > 0.0).mean(),
            "pz<0 frac": (df["pz"] < 0.0).mean(),
            "<theta_z> deg": df["theta_deg"].mean(),
            "median theta_z deg": df["theta_deg"].median(),
            "theta_z q10 deg": df["theta_deg"].quantile(0.10),
            "theta_z q90 deg": df["theta_deg"].quantile(0.90),
            "<phi> deg circular": circular_mean_deg(df["phi_deg"]),
            "|<exp(i phi)>|": circular_resultant(df["phi_deg"]),
            "<p_abs> beta_gamma": df["p_abs"].mean(),
            "median p_abs beta_gamma": df["p_abs"].median(),
            "<K> keV": df["kinetic_keV"].mean(),
            "median K keV": df["kinetic_keV"].median(),
        }
        for quantity, value in values.items():
            rows.append({"species": species, "quantity": quantity, "value": value})
    return pd.DataFrame(rows)


def summarize_pairs(electrons: pd.DataFrame, positrons: pd.DataFrame) -> dict[str, float]:
    """Compute event-wise electron/positron opening-angle diagnostics."""
    ue = electrons[["ux", "uy", "uz"]].to_numpy(float)
    up = positrons[["ux", "uy", "uz"]].to_numpy(float)
    pe = electrons[["px", "py", "pz"]].to_numpy(float)
    pp = positrons[["px", "py", "pz"]].to_numpy(float)

    dot = np.sum(ue * up, axis=1)
    opening = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))
    net = pe + pp
    net_norm = np.linalg.norm(net, axis=1)
    p_sum = electrons["p_abs"].to_numpy(float) + positrons["p_abs"].to_numpy(float)

    return {
        "<u_e dot u_p>": float(dot.mean()),
        "median opening angle deg": float(np.median(opening)),
        "opening q10 deg": float(np.quantile(opening, 0.10)),
        "opening q90 deg": float(np.quantile(opening, 0.90)),
        "fraction opening > 150 deg": float((opening > 150.0).mean()),
        "fraction opening < 30 deg": float((opening < 30.0).mean()),
        "<net px> beta_gamma": float(net[:, 0].mean()),
        "<net py> beta_gamma": float(net[:, 1].mean()),
        "<net pz> beta_gamma": float(net[:, 2].mean()),
        "median |p_e+p_p| beta_gamma": float(np.median(net_norm)),
        "median |p_e+p_p|/sum|p|": float(np.median(net_norm / p_sum)),
    }


def frames_have_time(frames: dict[str, pd.DataFrame]) -> bool:
    """Return true when each species dataframe has a creation-time column."""
    return all(TIME_COLUMN in df.columns for df in frames.values())


def summarize_time_species(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build per-species creation-time statistics in picoseconds."""
    rows = []
    for species, df in frames.items():
        t_ps = df[TIME_COLUMN].to_numpy(float) * 1.0e12
        values = {
            "N": len(t_ps),
            "t_min ps": np.min(t_ps),
            "t_max ps": np.max(t_ps),
            "t_span ps": np.ptp(t_ps),
            "t_mean ps": np.mean(t_ps),
            "t_std ps": np.std(t_ps, ddof=1),
        }
        for quantity, value in values.items():
            rows.append({"species": species, "quantity": quantity, "value": value})
    return pd.DataFrame(rows)


def summarize_time_pairs(electrons: pd.DataFrame, positrons: pd.DataFrame) -> dict[str, float]:
    """Compute pair-wise e-/e+ creation-time differences."""
    n_pairs = min(len(electrons), len(positrons))
    dt_s = (
        electrons[TIME_COLUMN].to_numpy(float)[:n_pairs]
        - positrons[TIME_COLUMN].to_numpy(float)[:n_pairs]
    )
    dt_ps = dt_s * 1.0e12
    abs_dt_ps = np.abs(dt_ps)
    return {
        "N pairs": float(n_pairs),
        "<t_e-t_p> ps": float(np.mean(dt_ps)),
        "std t_e-t_p ps": float(np.std(dt_ps, ddof=1)),
        "<|t_e-t_p|> ps": float(np.mean(abs_dt_ps)),
        "max |t_e-t_p| ps": float(np.max(abs_dt_ps)),
    }


def select_pair_index(electrons: pd.DataFrame, positrons: pd.DataFrame, requested: Optional[int]) -> int:
    """Return the requested pair index, or a presentation pair with a large theta."""
    n_pairs = min(len(electrons), len(positrons))
    if n_pairs == 0:
        raise ValueError("no e-/e+ pairs available")
    if requested is not None:
        if requested < 0 or requested >= n_pairs:
            raise ValueError(f"pair index {requested} outside [0, {n_pairs - 1}]")
        return requested

    ue = electrons[["ux", "uy", "uz"]].to_numpy(float)[:n_pairs]
    up = positrons[["ux", "uy", "uz"]].to_numpy(float)[:n_pairs]
    opening = np.degrees(np.arccos(np.clip(np.sum(ue * up, axis=1), -1.0, 1.0)))
    theta_e = electrons["theta_deg"].to_numpy(float)[:n_pairs]
    theta_p = positrons["theta_deg"].to_numpy(float)[:n_pairs]
    score = (
        np.abs(theta_e - 85.0)
        + 0.10 * np.abs(theta_p - 110.0)
        + 0.15 * np.abs(opening - 165.0)
    )
    return int(np.argmin(score))


def make_direction_figure(frames: dict[str, pd.DataFrame]):
    """Build the first-page momentum-direction diagnostic figure."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    colors = {"electron": "#2166ac", "positron": "#b2182b"}
    for species, df in frames.items():
        color = colors.get(species, None)
        axes[0, 0].hist(
            df["uz"], bins=60, histtype="step", lw=2, density=True, color=color, label=species
        )
        axes[0, 1].hist(
            df["theta_deg"],
            bins=60,
            histtype="step",
            lw=2,
            density=True,
            color=color,
            label=species,
        )
        axes[1, 0].hist(
            df["phi_deg"], bins=72, histtype="step", lw=2, density=True, color=color, label=species
        )
        axes[1, 1].scatter(df["ux"], df["uy"], s=4, alpha=0.25, color=color, label=species)

    axes[0, 0].set_xlabel(r"$u_z=p_z/|p|$")
    axes[0, 0].set_ylabel("density")
    axes[0, 0].set_title("Longitudinal momentum direction")
    axes[0, 1].set_xlabel(r"$\theta$ from +z [deg]")
    axes[0, 1].set_ylabel("density")
    axes[0, 1].set_title("Polar angle")
    axes[1, 0].set_xlabel(r"$\phi=\mathrm{atan2}(u_y,u_x)$ [deg]")
    axes[1, 0].set_ylabel("density")
    axes[1, 0].set_title("Azimuth")
    axes[1, 1].set_xlabel(r"$u_x$")
    axes[1, 1].set_ylabel(r"$u_y$")
    axes[1, 1].set_aspect("equal", adjustable="box")
    axes[1, 1].set_title("Transverse unit direction")
    for ax in axes.flat:
        ax.legend(frameon=False)
    return fig


def polar_arc(direction: np.ndarray, radius: float) -> np.ndarray:
    """Return arc points from +z to a unit direction at fixed radius."""
    transverse = np.array([direction[0], direction[1], 0.0])
    transverse_norm = np.linalg.norm(transverse)
    theta = float(np.arccos(np.clip(direction[2], -1.0, 1.0)))
    if transverse_norm < 1e-12 or theta < 1e-12:
        return np.empty((0, 3))
    rhat = transverse / transverse_norm
    zhat = np.array([0.0, 0.0, 1.0])
    samples = np.linspace(0.0, theta, 80)
    return radius * (np.sin(samples)[:, None] * rhat + np.cos(samples)[:, None] * zhat)


def azimuth_arc(phi: float, radius: float) -> np.ndarray:
    """Return xy-plane arc points from +x to the signed azimuth angle."""
    signed_phi = float(np.arctan2(np.sin(phi), np.cos(phi)))
    if abs(signed_phi) < 1e-12:
        return np.empty((0, 3))
    samples = np.linspace(0.0, signed_phi, 80)
    return radius * np.column_stack((np.cos(samples), np.sin(samples), np.zeros_like(samples)))


def opening_arc(a: np.ndarray, b: np.ndarray, radius: float) -> np.ndarray:
    """Return approximate great-circle arc points between two unit directions."""
    dot = float(np.clip(np.dot(a, b), -1.0, 1.0))
    alpha = float(np.arccos(dot))
    if alpha < 1e-12:
        return np.empty((0, 3))

    sin_alpha = np.sin(alpha)
    if abs(sin_alpha) < 1e-8:
        ortho = np.cross(a, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(ortho) < 1e-12:
            ortho = np.cross(a, np.array([0.0, 1.0, 0.0]))
        ortho /= np.linalg.norm(ortho)
        samples = np.linspace(0.0, alpha, 100)
        return radius * (np.cos(samples)[:, None] * a + np.sin(samples)[:, None] * ortho)

    samples = np.linspace(0.0, 1.0, 100)
    points = (
        np.sin((1.0 - samples) * alpha)[:, None] * a
        + np.sin(samples * alpha)[:, None] * b
    ) / sin_alpha
    return radius * points


def plot_arc(ax, points: np.ndarray, **kwargs) -> None:
    """Plot an arc when it contains visible points."""
    if len(points) > 0:
        ax.plot(points[:, 0], points[:, 1], points[:, 2], **kwargs)


def make_pair_geometry_figure(frames: dict[str, pd.DataFrame], pair_index: int):
    """Build the second-page coordinate-system and e-/e+ pair-angle figure."""
    import matplotlib.pyplot as plt

    electron = frames["electron"].iloc[pair_index]
    positron = frames["positron"].iloc[pair_index]
    ue = electron[["ux", "uy", "uz"]].to_numpy(float)
    up = positron[["ux", "uy", "uz"]].to_numpy(float)
    opening_deg = float(np.degrees(np.arccos(np.clip(np.dot(ue, up), -1.0, 1.0))))

    fig = plt.figure(figsize=(11, 8), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    axis_color = "#4d4d4d"
    for direction, label, offset in (
        (np.array([1.0, 0.0, 0.0]), "x", np.array([0.08, 0.0, 0.0])),
        (np.array([0.0, 1.0, 0.0]), "y", np.array([0.0, 0.08, 0.0])),
        (np.array([0.0, 0.0, 1.0]), "z", np.array([0.0, 0.0, 0.08])),
    ):
        ax.quiver(0, 0, 0, *direction, length=1.05, arrow_length_ratio=0.06, color=axis_color, lw=1.5)
        label_position = 1.08 * direction + offset
        ax.text(*label_position, label, color=axis_color, fontsize=14, weight="bold")

    vectors = [
        (ue, r"$e^-$", "#2166ac", electron),
        (up, r"$e^+$", "#b2182b", positron),
    ]
    for direction, label, color, row in vectors:
        ax.quiver(0, 0, 0, *direction, length=0.9, normalize=True, arrow_length_ratio=0.10, color=color, lw=3)
        ax.plot([0, 0.9 * direction[0]], [0, 0.9 * direction[1]], [0, 0], "--", color=color, alpha=0.5)
        ax.plot(
            [0.9 * direction[0], 0.9 * direction[0]],
            [0.9 * direction[1], 0.9 * direction[1]],
            [0.0, 0.9 * direction[2]],
            ":",
            color=color,
            alpha=0.5,
        )
        ax.text(*(0.98 * direction), label, color=color, fontsize=16, weight="bold")
        phi = np.radians(float(row["phi_deg"]))
        show_polar_arc = label != r"$e^+$"
        if show_polar_arc:
            plot_arc(ax, polar_arc(direction, 0.45), color=color, lw=2)
        plot_arc(ax, azimuth_arc(phi, 0.28), color=color, lw=2, ls="--")
        polar_label_position = polar_arc(direction, 0.52) if show_polar_arc else np.empty((0, 3))
        if len(polar_label_position) > 0:
            ax.text(*polar_label_position[len(polar_label_position) // 2], r"$\theta$", color=color, fontsize=12)
        az_label_position = azimuth_arc(phi, 0.33)
        if len(az_label_position) > 0:
            ax.text(*az_label_position[len(az_label_position) // 2], r"$\phi$", color=color, fontsize=12)
        energy_label_position = 0.78 * direction + np.array(
            [0.0, 0.0, 0.12 if direction[2] >= 0.0 else -0.16]
        )
        if label == r"$e^+$":
            energy_label_position += np.array([-0.08, -0.14, 0.28])
        energy_label_position = np.clip(energy_label_position, -0.82, 0.82)
        ax.text(
            *energy_label_position,
            f"K={row['kinetic_keV']:.1f} keV",
            color=color,
            fontsize=10,
            ha="center",
        )
    info = (
        f"Pair index: {pair_index}\n"
        f"opening alpha: {opening_deg:.2f} deg\n"
        f"theta e- / e+: {electron['theta_deg']:.2f} / {positron['theta_deg']:.2f} deg\n"
        f"phi e- / e+: {electron['phi_deg']:.2f} / {positron['phi_deg']:.2f} deg\n"
        f"|p| e- / e+: {electron['p_abs']:.4f} / {positron['p_abs']:.4f} beta*gamma\n"
        f"K e- / e+: {electron['kinetic_keV']:.1f} / {positron['kinetic_keV']:.1f} keV"
    )
    ax.text2D(
        0.03,
        0.97,
        info,
        transform=ax.transAxes,
        va="top",
        fontsize=11,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "#bdbdbd", "alpha": 0.92},
    )

    ax.set_title("Representative e-/e+ Pair Momentum Geometry", fontsize=16, pad=18)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_zlim(-1.05, 1.05)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=23, azim=38)
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["linewidth"] = 0
        axis._axinfo["grid"]["color"] = (1, 1, 1, 0)
    return fig


def make_time_figure(frames: dict[str, pd.DataFrame]):
    """Build creation-time histogram diagnostics for time-preserving FROMFILE data."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)
    colors = {"electron": "#2166ac", "positron": "#b2182b"}
    linestyles = {"electron": "-", "positron": "--"}
    for species, df in frames.items():
        color = colors.get(species, None)
        axes[0].hist(
            df[TIME_COLUMN] * 1.0e12,
            bins=60,
            histtype="step",
            lw=2.4,
            density=False,
            color=color,
            linestyle=linestyles.get(species, "-"),
            label=species,
        )

    dt_ps = (
        frames["electron"][TIME_COLUMN].to_numpy(float)
        - frames["positron"][TIME_COLUMN].to_numpy(float)
    ) * 1.0e12
    axes[1].hist(dt_ps, bins=60, histtype="stepfilled", color="#4d4d4d", alpha=0.75)
    axes[0].set_xlabel("creation time t [ps]")
    axes[0].set_ylabel("particles")
    axes[0].set_title("fort98 Creation-Time Distribution")
    axes[0].legend(frameon=False)
    axes[1].set_xlabel(r"$t_{e^-}-t_{e^+}$ [ps]")
    axes[1].set_ylabel("pairs")
    axes[1].set_title("Pair-Wise Creation-Time Difference")
    return fig


def write_plots(frames: dict[str, pd.DataFrame], out: Path, pair_index: Optional[int]) -> list[Path]:
    """Write PNG pages plus a multi-page PDF diagnostic."""
    cache_dir = out.parent / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 220,
            "font.size": 11,
            "axes.grid": True,
            "grid.alpha": 0.25,
        }
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    selected_pair = select_pair_index(frames["electron"], frames["positron"], pair_index)
    geometry_out = out.with_name("pair_momentum_pair_geometry.png")
    time_out = out.with_name("pair_momentum_time_histogram.png")
    pdf_out = out.with_suffix(".pdf")

    direction_fig = make_direction_figure(frames)
    geometry_fig = make_pair_geometry_figure(frames, selected_pair)
    time_fig = make_time_figure(frames) if frames_have_time(frames) else None

    direction_fig.savefig(out, bbox_inches="tight")
    geometry_fig.savefig(geometry_out, bbox_inches="tight")
    if time_fig is not None:
        time_fig.savefig(time_out, bbox_inches="tight")
    with PdfPages(pdf_out) as pdf:
        pdf.savefig(direction_fig, bbox_inches="tight")
        pdf.savefig(geometry_fig, bbox_inches="tight")
        if time_fig is not None:
            pdf.savefig(time_fig, bbox_inches="tight")

    plt.close(direction_fig)
    plt.close(geometry_fig)
    written = [out, geometry_out]
    if time_fig is not None:
        plt.close(time_fig)
        written.append(time_out)
    written.append(pdf_out)
    return written


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--electron-file", default="gamma_gamma_electrons.fromfile")
    parser.add_argument("--positron-file", default="gamma_gamma_positrons.fromfile")
    parser.add_argument("--output-dir", default="data")
    parser.add_argument("--plot", action="store_true", help="write PNG pages and a two-page PDF")
    parser.add_argument(
        "--pair-index",
        type=int,
        default=None,
        help="pair index for the geometry page; defaults to median opening-angle pair",
    )
    args = parser.parse_args()

    electron_file = Path(args.electron_file)
    positron_file = Path(args.positron_file)
    output_dir = Path(args.output_dir)

    frames = {
        "electron": add_momentum_columns(read_fromfile(electron_file)),
        "positron": add_momentum_columns(read_fromfile(positron_file)),
    }

    summary = summarize_species(frames)
    pivot = summary.pivot(index="quantity", columns="species", values="value")
    print(pivot.to_string(float_format=lambda x: f"{x:.6g}"))

    pair_summary = summarize_pairs(frames["electron"], frames["positron"])
    print("\nPair-wise e-/e+ direction relation:")
    for key, value in pair_summary.items():
        print(f"  {key:32s} = {value:.6g}")

    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "pair_momentum_species_summary.csv", index=False)
    pd.DataFrame([pair_summary]).to_csv(output_dir / "pair_momentum_pair_summary.csv", index=False)

    if frames_have_time(frames):
        time_summary = summarize_time_species(frames)
        time_pivot = time_summary.pivot(index="quantity", columns="species", values="value")
        print("\nCreation-time diagnostics:")
        print(time_pivot.to_string(float_format=lambda x: f"{x:.6g}"))

        time_pair_summary = summarize_time_pairs(frames["electron"], frames["positron"])
        print("\nPair-wise e-/e+ creation-time relation:")
        for key, value in time_pair_summary.items():
            print(f"  {key:32s} = {value:.6g}")

        time_summary.to_csv(output_dir / "pair_time_species_summary.csv", index=False)
        pd.DataFrame([time_pair_summary]).to_csv(output_dir / "pair_time_pair_summary.csv", index=False)

    if args.plot:
        out = output_dir / "pair_momentum_directions.png"
        written = write_plots(frames, out, args.pair_index)
        print("\nWrote")
        for path in written:
            print(f"  {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
