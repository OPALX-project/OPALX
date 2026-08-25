#!/usr/bin/env python3
"""Analytic boosted-Gaussian witness model matching the current OPALX setup.

This script uses the triaxial rigid Gaussian field evaluator from
``sandbox/track12particles/track12particles.py`` and tracks the six e- plus six
e+ FROMFILE witnesses through the same timing window used by the clean OPALX
BeamBeam run.
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

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_OPALX_RUN = (
    ROOT
    / "sandbox"
    / "Drift-Experiment"
    / "one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs_extended_bb_nobound"
)

ELECTRON_COLOR = "#c00000"
POSITRON_COLOR = "#0017b8"
LINE_STYLES = {
    1: "-",
    2: (0, (8, 5)),
    3: (0, (5, 5)),
    4: (0, (9, 6)),
    5: (0, (2, 3)),
    6: (0, (10, 4, 2, 4)),
}


@dataclass(frozen=True)
class Setup:
    elementary_charge_C: float
    electron_rest_eV: float
    c_light_m_per_s: float
    primary_kinetic_GeV: float
    primary_gamma: float
    primary_beta: float
    primary_electrons_per_bunch: float
    primary_charge_C: float
    primary_macroparticles: int
    sigma_x_m: float
    sigma_y_m: float
    sigma_z_m: float
    witness_particles_per_species: int
    witness_bunch_charge_C: float
    witness_t0_s: float
    witness_r0z_m: float
    bb_start_s_m: float
    bb_ip_s_m: float
    bb_length_m: float
    bb_end_s_m: float
    primary_edge_sigma: float
    track12_first_ct_m: float
    track12_edge_reference_s_m: float
    primary_centroid_at_witness_t0_m: float
    primary_source_r0z_m: float
    near_ip_active_time_s: float
    post_retire_observation_time_s: float
    primary_retire_time_s: float
    fine_dt_s: float
    ps_dump_freq: int
    aperture_radius_m: float


def default_setup(source_charge_scale: float) -> Setup:
    elementary_charge = 1.602_176_634e-19
    electron_rest_eV = 510_998.95
    c_light = 299_792_458.0
    emass_GeV = electron_rest_eV * 1.0e-9

    primary_kinetic_GeV = 0.245
    primary_gamma = (primary_kinetic_GeV + emass_GeV) / emass_GeV
    primary_beta = math.sqrt(1.0 - 1.0 / (primary_gamma * primary_gamma))

    primary_electrons = 1.25e10 * source_charge_scale
    primary_charge_C = -primary_electrons * elementary_charge
    sigma_xy = 1.944_325_075_701e-6
    sigma_z = 0.6e-3

    witness_t0 = 4.0e-12
    bb_start_s = 1.0e-3
    bb_ip_s = 6.0e-3
    bb_length = 1.9e-2
    primary_edge_sigma = 3.0
    track12_first_ct = -1.5 * sigma_z
    track12_edge_reference_s = bb_ip_s + track12_first_ct
    primary_centroid_at_witness_t0 = track12_edge_reference_s - primary_edge_sigma * sigma_z
    primary_source_r0z = primary_centroid_at_witness_t0 - primary_beta * c_light * witness_t0
    near_ip_active_time = 5.0e-11
    post_retire_observation_time = 1.0e-12

    witness_particles_per_species = 6
    witness_bunch_charge = witness_particles_per_species * elementary_charge

    return Setup(
        elementary_charge_C=elementary_charge,
        electron_rest_eV=electron_rest_eV,
        c_light_m_per_s=c_light,
        primary_kinetic_GeV=primary_kinetic_GeV,
        primary_gamma=primary_gamma,
        primary_beta=primary_beta,
        primary_electrons_per_bunch=primary_electrons,
        primary_charge_C=primary_charge_C,
        primary_macroparticles=400_000,
        sigma_x_m=sigma_xy,
        sigma_y_m=sigma_xy,
        sigma_z_m=sigma_z,
        witness_particles_per_species=witness_particles_per_species,
        witness_bunch_charge_C=witness_bunch_charge,
        witness_t0_s=witness_t0,
        witness_r0z_m=bb_ip_s,
        bb_start_s_m=bb_start_s,
        bb_ip_s_m=bb_ip_s,
        bb_length_m=bb_length,
        bb_end_s_m=bb_start_s + bb_length,
        primary_edge_sigma=primary_edge_sigma,
        track12_first_ct_m=track12_first_ct,
        track12_edge_reference_s_m=track12_edge_reference_s,
        primary_centroid_at_witness_t0_m=primary_centroid_at_witness_t0,
        primary_source_r0z_m=primary_source_r0z,
        near_ip_active_time_s=near_ip_active_time,
        post_retire_observation_time_s=post_retire_observation_time,
        primary_retire_time_s=witness_t0 + near_ip_active_time,
        fine_dt_s=5.0e-15,
        ps_dump_freq=20,
        aperture_radius_m=5.0e-4,
    )


def configure_matplotlib(output_dir: Path) -> None:
    cache_dir = output_dir / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sorted_steps(h5file: h5py.File) -> list[str]:
    return sorted(h5file.keys(), key=lambda key: int(key.split("#")[1]))


def reference_z(group: h5py.Group) -> float:
    if "RefPartR" in group.attrs:
        return float(group.attrs["RefPartR"][2])
    if "SPOS" in group.attrs:
        return float(group.attrs["SPOS"][0])
    return 0.0


def group_time_s(group: h5py.Group) -> float:
    if "TIME" in group.attrs:
        return float(group.attrs["TIME"][0])
    if "T" in group.attrs:
        return float(group.attrs["T"][0])
    raise KeyError("H5 group does not contain TIME or T")


def h5_step_index(step_name: str) -> int:
    return int(step_name.split("#", maxsplit=1)[1])


def read_fromfile(
    path: Path,
    *,
    species: str,
    kind: int,
    charge_units: int,
    setup: Setup,
) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as stream:
        expected_n = int(stream.readline().strip())
        columns = stream.readline().strip().split()
        rows = pd.read_csv(stream, sep=r"\s+", names=columns)

    if len(rows) != expected_n:
        raise ValueError(f"{path}: header says {expected_n} particles, read {len(rows)}")
    if len(rows) != setup.witness_particles_per_species:
        raise ValueError(
            f"{path}: expected {setup.witness_particles_per_species} particles for the OPALX setup, read {len(rows)}"
        )

    rows = rows.copy()
    rows["pair"] = np.arange(1, len(rows) + 1, dtype=int)
    rows["kind"] = kind
    rows["species"] = species
    rows["charge_units"] = charge_units
    rows["name"] = [f"Te{'-' if kind == 2 else '+'}{pair}" for pair in rows["pair"]]
    rows["x_m"] = rows["x"]
    rows["y_m"] = rows["y"]
    rows["z_abs_m"] = setup.witness_r0z_m + rows["z"]
    rows["s_minus_ip_m"] = rows["z_abs_m"] - setup.bb_ip_s_m
    rows["px"] = rows["px"]
    rows["py"] = rows["py"]
    rows["pz"] = rows["pz"]
    rows["witness_bunch_charge_C"] = math.copysign(setup.witness_bunch_charge_C, charge_units)
    return rows


def load_initial_witnesses(input_dir: Path, setup: Setup) -> pd.DataFrame:
    electrons = read_fromfile(
        input_dir / "track12_electrons.fromfile",
        species="electron",
        kind=2,
        charge_units=-1,
        setup=setup,
    )
    positrons = read_fromfile(
        input_dir / "track12_positrons.fromfile",
        species="positron",
        kind=3,
        charge_units=+1,
        setup=setup,
    )
    return pd.concat([electrons, positrons], ignore_index=True)


def make_source(track12, setup: Setup):
    return track12.RigidAnisotropicGaussianSource(
        "c0",
        charge_C=setup.primary_charge_C,
        sigma_lab_m=(setup.sigma_x_m, setup.sigma_y_m, setup.sigma_z_m),
        beta_z=setup.primary_beta,
        center_t0_m=np.array([0.0, 0.0, setup.primary_source_r0z_m], dtype=float),
        t0_s=0.0,
    )


def source_is_active(source, time_s: float, setup: Setup) -> bool:
    center_z = float(source.center(time_s)[2])
    in_beambeam = setup.bb_start_s_m <= center_z <= setup.bb_end_s_m
    before_retire = time_s <= setup.primary_retire_time_s + 1.0e-21
    return bool(in_beambeam and before_retire)


def make_field_function(track12, setup: Setup):
    def field_function(position_m: np.ndarray, time_s: float, sources) -> tuple[np.ndarray, np.ndarray]:
        active_sources = tuple(source for source in sources if source_is_active(source, time_s, setup))
        if not active_sources:
            return np.zeros(3), np.zeros(3)
        return track12.anisotropic_total_lab_fields(position_m, time_s, active_sources)

    return field_function


def load_opalx_target_times(run_dir: Path, setup: Setup) -> list[float]:
    path = run_dir / "spacecharge_drift_withness_c1.h5"
    if not path.exists():
        final_time = setup.primary_retire_time_s + setup.post_retire_observation_time_s + 5.0e-14
        cadence = setup.ps_dump_freq * setup.fine_dt_s
        generated = np.arange(setup.witness_t0_s + cadence, final_time + 0.5 * cadence, cadence)
        return [float(t) for t in generated]
    with h5py.File(path, "r") as h5:
        return [group_time_s(h5[key]) for key in sorted_steps(h5)]


def record_particle(
    row: dict[str, object],
    *,
    position: np.ndarray,
    momentum: np.ndarray,
    e_field: np.ndarray,
    b_field: np.ndarray,
    time_s: float,
    local_step: int,
    source,
    setup: Setup,
) -> dict[str, object]:
    gamma = math.sqrt(1.0 + float(np.dot(momentum, momentum)))
    source_center = source.center(time_s)
    return {
        "model": "analytic Gaussian c0",
        "name": row["name"],
        "pair": int(row["pair"]),
        "kind": int(row["kind"]),
        "species": row["species"],
        "charge_units": int(row["charge_units"]),
        "witness_bunch_charge_C": float(row["witness_bunch_charge_C"]),
        "step_since_t0": local_step,
        "time_s": time_s,
        "time_fs": int(round(time_s * 1.0e15)),
        "time_ps": time_s * 1.0e12,
        "time_since_t0_fs": (time_s - setup.witness_t0_s) * 1.0e15,
        "time_since_t0_ps": (time_s - setup.witness_t0_s) * 1.0e12,
        "source_active": source_is_active(source, time_s, setup),
        "source_z_m": float(source_center[2]),
        "source_plus_3sigma_z_m": float(source_center[2] + 3.0 * setup.sigma_z_m),
        "source_charge_C": setup.primary_charge_C,
        "x_m": float(position[0]),
        "y_m": float(position[1]),
        "z_abs_m": float(position[2]),
        "s_minus_ip_m": float(position[2] - setup.bb_ip_s_m),
        "x_um": float(position[0] * 1.0e6),
        "y_um": float(position[1] * 1.0e6),
        "s_minus_ip_mm": float((position[2] - setup.bb_ip_s_m) * 1.0e3),
        "px": float(momentum[0]),
        "py": float(momentum[1]),
        "pz": float(momentum[2]),
        "gamma": gamma,
        "kinetic_keV": (gamma - 1.0) * setup.electron_rest_eV * 1.0e-3,
        "Ex_V_per_m": float(e_field[0]),
        "Ey_V_per_m": float(e_field[1]),
        "Ez_V_per_m": float(e_field[2]),
        "Bx_T": float(b_field[0]),
        "By_T": float(b_field[1]),
        "Bz_T": float(b_field[2]),
        "E_abs_V_per_m": float(np.linalg.norm(e_field)),
        "B_abs_T": float(np.linalg.norm(b_field)),
    }


def track_analytic(
    initial: pd.DataFrame,
    target_times_s: list[float],
    track12,
    source,
    setup: Setup,
) -> pd.DataFrame:
    target_steps = {0}
    for time_s in target_times_s:
        target_steps.add(int(round((time_s - setup.witness_t0_s) / setup.fine_dt_s)))
    if min(target_steps) < 0:
        raise ValueError("target OPALX time precedes witness T0")

    max_step = max(target_steps)
    sources = (source,)
    field_function = make_field_function(track12, setup)
    particles = []
    for row in initial.to_dict(orient="records"):
        particles.append(
            {
                "meta": row,
                "position": np.array([row["x_m"], row["y_m"], row["z_abs_m"]], dtype=float),
                "momentum": np.array([row["px"], row["py"], row["pz"]], dtype=float),
            }
        )

    records: list[dict[str, object]] = []
    for local_step in range(max_step + 1):
        time_s = setup.witness_t0_s + local_step * setup.fine_dt_s
        if local_step in target_steps:
            for particle in particles:
                e_field, b_field = field_function(particle["position"], time_s, sources)
                records.append(
                    record_particle(
                        particle["meta"],
                        position=particle["position"],
                        momentum=particle["momentum"],
                        e_field=e_field,
                        b_field=b_field,
                        time_s=time_s,
                        local_step=local_step,
                        source=source,
                        setup=setup,
                    )
                )

        if local_step == max_step:
            break

        for particle in particles:
            position, momentum, _e_field, _b_field = track12.advance_with_fields(
                particle["position"],
                particle["momentum"],
                float(particle["meta"]["charge_units"]),
                time_s,
                setup.fine_dt_s,
                field_function,
                sources,
            )
            particle["position"] = position
            particle["momentum"] = momentum

    return pd.DataFrame(records).sort_values(["species", "pair", "time_s"]).reset_index(drop=True)


def load_opalx_witness_h5(run_dir: Path, setup: Setup) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    files = [
        (run_dir / "spacecharge_drift_withness_c1.h5", "electron", 2, -1),
        (run_dir / "spacecharge_drift_withness_c2.h5", "positron", 3, +1),
    ]
    for path, species, kind, charge_units in files:
        if not path.exists():
            continue
        with h5py.File(path, "r") as h5:
            for step_name in sorted_steps(h5):
                group = h5[step_name]
                time_s = group_time_s(group)
                ref_z = reference_z(group)
                ids = group["id"][:] if "id" in group else np.arange(len(group["x"]))
                for index, particle_id in enumerate(ids):
                    pair = int(particle_id) + 1
                    x = float(group["x"][index])
                    y = float(group["y"][index])
                    z_abs = ref_z + float(group["z"][index])
                    momentum = np.array(
                        [float(group["px"][index]), float(group["py"][index]), float(group["pz"][index])],
                        dtype=float,
                    )
                    gamma = math.sqrt(1.0 + float(np.dot(momentum, momentum)))
                    e_field = np.array(
                        [float(group["Ex"][index]), float(group["Ey"][index]), float(group["Ez"][index])],
                        dtype=float,
                    )
                    b_field = np.array(
                        [float(group["Bx"][index]), float(group["By"][index]), float(group["Bz"][index])],
                        dtype=float,
                    )
                    rows.append(
                        {
                            "model": "OPALX H5 BeamBeam",
                            "h5_file": str(path),
                            "h5_step": h5_step_index(step_name),
                            "global_step": int(group.attrs["GlobalTrackStep"][0])
                            if "GlobalTrackStep" in group.attrs
                            else np.nan,
                            "pair": pair,
                            "kind": kind,
                            "species": species,
                            "charge_units": charge_units,
                            "witness_bunch_charge_C": float(group.attrs["CHARGE"][0])
                            if "CHARGE" in group.attrs
                            else math.copysign(setup.witness_bunch_charge_C, charge_units),
                            "time_s": time_s,
                            "time_fs": int(round(time_s * 1.0e15)),
                            "time_ps": time_s * 1.0e12,
                            "time_since_t0_fs": (time_s - setup.witness_t0_s) * 1.0e15,
                            "time_since_t0_ps": (time_s - setup.witness_t0_s) * 1.0e12,
                            "x_m": x,
                            "y_m": y,
                            "z_abs_m": z_abs,
                            "s_minus_ip_m": z_abs - setup.bb_ip_s_m,
                            "x_um": x * 1.0e6,
                            "y_um": y * 1.0e6,
                            "s_minus_ip_mm": (z_abs - setup.bb_ip_s_m) * 1.0e3,
                            "px": float(momentum[0]),
                            "py": float(momentum[1]),
                            "pz": float(momentum[2]),
                            "gamma": gamma,
                            "kinetic_keV": (gamma - 1.0) * setup.electron_rest_eV * 1.0e-3,
                            "Ex_V_per_m": float(e_field[0]),
                            "Ey_V_per_m": float(e_field[1]),
                            "Ez_V_per_m": float(e_field[2]),
                            "Bx_T": float(b_field[0]),
                            "By_T": float(b_field[1]),
                            "Bz_T": float(b_field[2]),
                            "E_abs_V_per_m": float(np.linalg.norm(e_field)),
                            "B_abs_T": float(np.linalg.norm(b_field)),
                        }
                    )
    return pd.DataFrame(rows).sort_values(["species", "pair", "time_s"]).reset_index(drop=True)


def make_comparison(analytic: pd.DataFrame, opalx: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if opalx.empty:
        return pd.DataFrame(), pd.DataFrame()

    analytic_samples = analytic.loc[analytic["time_fs"].isin(opalx["time_fs"].unique())].copy()
    merged = analytic_samples.merge(
        opalx,
        on=["species", "pair", "kind", "time_fs"],
        suffixes=("_analytic", "_opalx"),
    )
    difference_columns = [
        "x_um",
        "y_um",
        "s_minus_ip_mm",
        "px",
        "py",
        "pz",
        "kinetic_keV",
        "Ex_V_per_m",
        "Ey_V_per_m",
        "Ez_V_per_m",
        "E_abs_V_per_m",
    ]
    for column in difference_columns:
        merged[f"d_{column}"] = merged[f"{column}_analytic"] - merged[f"{column}_opalx"]

    summary_rows = []
    for column in difference_columns:
        diff = merged[f"d_{column}"].to_numpy(dtype=float)
        opalx_values = merged[f"{column}_opalx"].to_numpy(dtype=float)
        denom = np.linalg.norm(opalx_values)
        summary_rows.append(
            {
                "quantity": column,
                "samples": len(diff),
                "mean_difference": float(np.mean(diff)),
                "rms_difference": float(np.sqrt(np.mean(diff * diff))),
                "max_abs_difference": float(np.max(np.abs(diff))),
                "relative_l2_difference": float(np.linalg.norm(diff) / denom) if denom else np.nan,
            }
        )
    summary = pd.DataFrame(summary_rows)
    return merged, summary


def plot_model_groups(ax, frame: pd.DataFrame, *, x_column: str, y_column: str, alpha: float = 1.0) -> None:
    for (species, pair), group in frame.groupby(["species", "pair"], sort=False):
        color = ELECTRON_COLOR if species == "electron" else POSITRON_COLOR
        group = group.sort_values("time_s")
        ax.plot(
            group[x_column],
            group[y_column],
            color=color,
            ls=LINE_STYLES.get(int(pair), "-"),
            lw=1.25,
            alpha=alpha,
        )


def plot_side_by_side(analytic: pd.DataFrame, opalx: pd.DataFrame, output: Path, *, x_axis: str) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    if x_axis == "s":
        x_column = "s_minus_ip_mm"
        xlabel = "S - S_IP [mm]"
    elif x_axis == "t":
        x_column = "time_since_t0_ps"
        xlabel = "t - T0 [ps]"
    else:
        raise ValueError(f"unsupported x axis {x_axis!r}")

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2), sharey=True, constrained_layout=True)
    if not opalx.empty:
        plot_model_groups(axes[0], opalx, x_column=x_column, y_column="x_um")
    plot_model_groups(axes[1], analytic, x_column=x_column, y_column="x_um")

    axes[0].set_title("OPALX H5 BeamBeam")
    axes[1].set_title("analytic Gaussian c0")
    for ax in axes:
        ax.axhline(0.0, color="black", lw=0.9, ls=(0, (2, 6)))
        ax.set_xlabel(xlabel)
        ax.grid(True, alpha=0.25)
        ax.minorticks_on()
        ax.tick_params(direction="out", top=True, right=True, which="both")
    axes[0].set_ylabel("x [um]")
    fig.suptitle("track12 witnesses: x trajectory")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_overlay(analytic: pd.DataFrame, opalx: pd.DataFrame, output: Path, *, x_axis: str) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    if x_axis == "s":
        x_column = "s_minus_ip_mm"
        xlabel = "S - S_IP [mm]"
    elif x_axis == "t":
        x_column = "time_since_t0_ps"
        xlabel = "t - T0 [ps]"
    else:
        raise ValueError(f"unsupported x axis {x_axis!r}")

    fig, ax = plt.subplots(figsize=(8.6, 5.4), constrained_layout=True)
    for frame, model_style, model_lw in ((opalx, "-", 1.4), (analytic, "--", 1.2)):
        if frame.empty:
            continue
        for (species, pair), group in frame.groupby(["species", "pair"], sort=False):
            color = ELECTRON_COLOR if species == "electron" else POSITRON_COLOR
            group = group.sort_values("time_s")
            ax.plot(
                group[x_column],
                group["x_um"],
                color=color,
                ls=model_style,
                lw=model_lw,
                alpha=0.92,
            )

    ax.axhline(0.0, color="black", lw=0.9, ls=(0, (2, 6)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("x [um]")
    ax.set_title("analytic Gaussian c0 vs OPALX")
    ax.grid(True, alpha=0.25)
    ax.minorticks_on()
    species_handles = [
        plt.Line2D([0], [0], color=ELECTRON_COLOR, lw=2.0, label="e-"),
        plt.Line2D([0], [0], color=POSITRON_COLOR, lw=2.0, label="e+"),
    ]
    model_handles = [
        plt.Line2D([0], [0], color="0.2", lw=1.5, ls="-", label="OPALX"),
        plt.Line2D([0], [0], color="0.2", lw=1.5, ls="--", label="analytic"),
    ]
    first = ax.legend(handles=species_handles, loc="upper left", fontsize=9, title="species")
    ax.add_artist(first)
    ax.legend(handles=model_handles, loc="upper right", fontsize=9, title="model")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_field_timing(analytic: pd.DataFrame, opalx: pd.DataFrame, output: Path, setup: Setup) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    frames = []
    if not opalx.empty:
        frames.append(opalx.assign(model_short="OPALX"))
    frames.append(analytic.assign(model_short="analytic"))
    combined = pd.concat(frames, ignore_index=True)

    grouped = (
        combined.groupby(["model_short", "time_fs"], as_index=False)
        .agg(
            time_since_t0_ps=("time_since_t0_ps", "first"),
            E_max_V_per_m=("E_abs_V_per_m", "max"),
            E_median_V_per_m=("E_abs_V_per_m", "median"),
        )
        .sort_values(["model_short", "time_fs"])
    )

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True, constrained_layout=True)
    for model, group in grouped.groupby("model_short", sort=False):
        style = "-" if model == "OPALX" else "--"
        color = "tab:blue" if model == "OPALX" else "tab:orange"
        axes[0].semilogy(
            group["time_since_t0_ps"],
            group["E_max_V_per_m"],
            color=color,
            ls=style,
            lw=1.6,
            label=f"{model} max |E|",
        )
        axes[0].semilogy(
            group["time_since_t0_ps"],
            group["E_median_V_per_m"],
            color=color,
            ls=(0, (2, 3)) if model == "OPALX" else (0, (5, 4)),
            lw=1.1,
            alpha=0.85,
            label=f"{model} median |E|",
        )

    time_grid = np.linspace(setup.witness_t0_s, setup.primary_retire_time_s, 300)
    source_z_mm = (
        setup.primary_source_r0z_m
        + setup.primary_beta * setup.c_light_m_per_s * time_grid
        - setup.bb_ip_s_m
    ) * 1.0e3
    source_edge_mm = source_z_mm + 3.0 * setup.sigma_z_m * 1.0e3
    axes[1].plot((time_grid - setup.witness_t0_s) * 1.0e12, source_z_mm, color="tab:green", lw=1.5, label="c0 centroid")
    axes[1].plot(
        (time_grid - setup.witness_t0_s) * 1.0e12,
        source_edge_mm,
        color="tab:green",
        lw=1.1,
        ls="--",
        label="c0 +3 sigma_z edge",
    )

    for ax in axes:
        ax.axvline(
            (setup.primary_retire_time_s - setup.witness_t0_s) * 1.0e12,
            color="black",
            ls="--",
            lw=1.0,
            label="c0 retire" if ax is axes[0] else None,
        )
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(loc="best", fontsize=8)
    axes[0].set_ylabel("|E| [V/m]")
    axes[0].set_title("sampled witness field timing")
    axes[1].set_xlabel("t - T0 [ps]")
    axes[1].set_ylabel("source S - S_IP [mm]")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_differences(comparison: pd.DataFrame, output: Path) -> None:
    if comparison.empty:
        return
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.2), sharex=True, constrained_layout=True)
    quantities = [
        ("d_x_um", "analytic - OPALX x [um]"),
        ("d_y_um", "analytic - OPALX y [um]"),
        ("d_E_abs_V_per_m", "analytic - OPALX |E| [V/m]"),
    ]
    for ax, (column, ylabel) in zip(axes, quantities, strict=True):
        for (species, pair), group in comparison.groupby(["species", "pair"], sort=False):
            color = ELECTRON_COLOR if species == "electron" else POSITRON_COLOR
            group = group.sort_values("time_s_analytic")
            ax.plot(
                group["time_since_t0_ps_analytic"],
                group[column],
                color=color,
                ls=LINE_STYLES.get(int(pair), "-"),
                lw=1.0,
                alpha=0.9,
            )
        ax.axhline(0.0, color="black", lw=0.8)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("t - T0 [ps]")
    species_handles = [
        plt.Line2D([0], [0], color=ELECTRON_COLOR, lw=2.0, label="e-"),
        plt.Line2D([0], [0], color=POSITRON_COLOR, lw=2.0, label="e+"),
    ]
    axes[0].legend(handles=species_handles, loc="best", fontsize=8)
    fig.suptitle("analytic minus OPALX on common H5 samples")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_outputs(
    *,
    output_dir: Path,
    setup: Setup,
    initial: pd.DataFrame,
    analytic: pd.DataFrame,
    opalx: pd.DataFrame,
    comparison: pd.DataFrame,
    comparison_summary: pd.DataFrame,
) -> list[Path]:
    paths = []
    setup_path = output_dir / "setup_summary.json"
    setup_path.write_text(json.dumps(asdict(setup), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    paths.append(setup_path)

    initial_path = output_dir / "witness_initial_conditions.csv"
    initial.to_csv(initial_path, index=False)
    paths.append(initial_path)

    analytic_path = output_dir / "analytic_witness_trajectory.csv"
    analytic.to_csv(analytic_path, index=False)
    paths.append(analytic_path)

    if not opalx.empty:
        opalx_path = output_dir / "opalx_witness_trajectory_sampled.csv"
        opalx.to_csv(opalx_path, index=False)
        paths.append(opalx_path)

    if not comparison.empty:
        comparison_path = output_dir / "analytic_vs_opalx_samples.csv"
        summary_path = output_dir / "analytic_vs_opalx_summary.csv"
        comparison.to_csv(comparison_path, index=False)
        comparison_summary.to_csv(summary_path, index=False)
        paths.extend([comparison_path, summary_path])

    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory containing track12_electrons.fromfile and track12_positrons.fromfile.",
    )
    parser.add_argument("--opalx-run-dir", type=Path, default=DEFAULT_OPALX_RUN)
    parser.add_argument(
        "--source-charge-scale",
        type=float,
        default=1.0,
        help="Scale applied to the c0 electron count. Negative values flip the analytic field sign.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    setup = default_setup(args.source_charge_scale)
    track12 = load_track12_module()
    source = make_source(track12, setup)
    initial = load_initial_witnesses(args.input_dir.resolve(), setup)
    target_times_s = load_opalx_target_times(args.opalx_run_dir, setup)

    analytic = track_analytic(initial, target_times_s, track12, source, setup)
    opalx = load_opalx_witness_h5(args.opalx_run_dir, setup)
    comparison, comparison_summary = make_comparison(analytic, opalx)

    written = write_outputs(
        output_dir=output_dir,
        setup=setup,
        initial=initial,
        analytic=analytic,
        opalx=opalx,
        comparison=comparison,
        comparison_summary=comparison_summary,
    )

    plot_paths = [
        output_dir / "analytic_vs_opalx_side_by_side_x_vs_s.png",
        output_dir / "analytic_vs_opalx_side_by_side_x_vs_t.png",
        output_dir / "analytic_vs_opalx_overlay_x_vs_s.png",
        output_dir / "analytic_vs_opalx_overlay_x_vs_t.png",
        output_dir / "analytic_field_timing.png",
    ]
    plot_side_by_side(analytic, opalx, plot_paths[0], x_axis="s")
    plot_side_by_side(analytic, opalx, plot_paths[1], x_axis="t")
    plot_overlay(analytic, opalx, plot_paths[2], x_axis="s")
    plot_overlay(analytic, opalx, plot_paths[3], x_axis="t")
    plot_field_timing(analytic, opalx, plot_paths[4], setup)
    if not comparison.empty:
        difference_plot = output_dir / "analytic_minus_opalx_differences.png"
        plot_differences(comparison, difference_plot)
        plot_paths.append(difference_plot)

    print(f"Output directory: {output_dir}")
    print(f"Analytic trajectory rows: {len(analytic)}")
    if not opalx.empty:
        print(f"OPALX trajectory rows: {len(opalx)}")
        print(f"Common comparison rows: {len(comparison)}")
    print("\nKey setup:")
    print(f"  c0 charge: {setup.primary_charge_C:.12e} C")
    print(f"  c0 beta: {setup.primary_beta:.12f}")
    print(f"  c0 sigma xy/z: {setup.sigma_x_m:.12e} m / {setup.sigma_z_m:.12e} m")
    print(f"  witness T0: {setup.witness_t0_s * 1.0e12:.6g} ps")
    print(f"  retire time: {setup.primary_retire_time_s * 1.0e12:.6g} ps")
    print(f"  c0 center at T0: {setup.primary_centroid_at_witness_t0_m * 1.0e3:.12g} mm")
    print(f"  c0 +3sigma edge at T0: {(setup.primary_centroid_at_witness_t0_m + 3.0 * setup.sigma_z_m) * 1.0e3:.12g} mm")
    print(f"  first witness S: {setup.track12_edge_reference_s_m * 1.0e3:.12g} mm")

    if not comparison_summary.empty:
        print("\nLargest comparison differences:")
        ordered = comparison_summary.sort_values("relative_l2_difference", ascending=False)
        print(ordered.head(8).to_string(index=False))

    print("\nWrote data:")
    for path in written:
        print(f"  {path}")
    print("Wrote plots:")
    for path in plot_paths:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
