#!/usr/bin/env python3
"""Manufactured Lorentz-transformed Gaussian source fields for one witness.

This script is intentionally independent of OPALX internals.  It evaluates the
closed-form rest-frame electric field of one or two spherical Gaussian charge
clouds, Lorentz-transforms the fields to the lab/reference frame, and applies
one OPALX-style Boris kick to a witness electron.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


C_LIGHT = 299_792_458.0
EPSILON_0 = 8.854_187_812_8e-12
ELEMENTARY_CHARGE_C = 1.602_176_634e-19
ELECTRON_REST_EV = 510_998.95
ELECTRON_CHARGE_UNITS = -1.0
DEFAULT_SOURCE_ELECTRONS = 1.25e10
DEFAULT_SOURCE_CHARGE_C = -DEFAULT_SOURCE_ELECTRONS * ELEMENTARY_CHARGE_C
SANDBOX_DIR = Path(__file__).resolve().parents[1]
NOTE_DIR = SANDBOX_DIR / "note"
DATA_DIR = SANDBOX_DIR / "data"


@dataclass(frozen=True)
class Source:
    """Uniformly moving spherical Gaussian source bunch."""

    name: str
    charge_C: float
    sigma_m: float
    beta: np.ndarray
    center_t0_m: np.ndarray
    t0_s: float = 0.0

    @property
    def gamma(self) -> float:
        beta2 = float(np.dot(self.beta, self.beta))
        if beta2 >= 1.0:
            raise ValueError(f"{self.name}: |beta| must be < 1")
        return 1.0 / math.sqrt(1.0 - beta2)

    @property
    def beta_hat(self) -> np.ndarray:
        beta_norm = float(np.linalg.norm(self.beta))
        if beta_norm == 0.0:
            return np.zeros(3)
        return self.beta / beta_norm

    def center(self, time_s: float) -> np.ndarray:
        """Return the lab-frame source center at time_s."""
        return self.center_t0_m + self.beta * C_LIGHT * (time_s - self.t0_s)


@dataclass(frozen=True)
class WitnessCase:
    """One manufactured witness-kick case."""

    name: str
    description: str
    position_m: np.ndarray
    momentum: np.ndarray
    sources: tuple[Source, ...]


@dataclass(frozen=True)
class WitnessParticle:
    """One witness particle integrated through the manufactured fields."""

    name: str
    charge_units: float
    position_m: np.ndarray
    momentum: np.ndarray


def beta_from_kinetic_energy(kinetic_MeV: float, rest_eV: float = ELECTRON_REST_EV) -> float:
    """Return beta for an electron with the given kinetic energy."""
    gamma = 1.0 + kinetic_MeV * 1.0e6 / rest_eV
    return math.sqrt(1.0 - 1.0 / (gamma * gamma))


def momentum_magnitude_from_kinetic_energy_keV(
    kinetic_keV: float, rest_eV: float = ELECTRON_REST_EV
) -> float:
    """Return |P|=beta*gamma from kinetic energy in keV."""
    gamma = 1.0 + kinetic_keV * 1.0e3 / rest_eV
    return math.sqrt(max(gamma * gamma - 1.0, 0.0))


def gamma_from_momentum(momentum: np.ndarray) -> float:
    """Return relativistic gamma from normalized momentum beta*gamma."""
    return math.sqrt(1.0 + float(np.dot(momentum, momentum)))


def velocity_from_momentum(momentum: np.ndarray) -> np.ndarray:
    """Return lab velocity from normalized momentum beta*gamma."""
    return C_LIGHT * momentum / gamma_from_momentum(momentum)


def split_parallel_perpendicular(vector: np.ndarray, direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split vector into components parallel/perpendicular to unit direction."""
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        return np.zeros(3), vector.copy()
    unit = direction / norm
    parallel = float(np.dot(vector, unit)) * unit
    return parallel, vector - parallel


def gaussian_rest_field(r_prime_m: np.ndarray, charge_C: float, sigma_m: float) -> tuple[float, np.ndarray]:
    """Return rest-frame potential and E field for a spherical Gaussian charge."""
    radius2 = float(np.dot(r_prime_m, r_prime_m))
    radius = math.sqrt(radius2)
    sigma2 = sigma_m * sigma_m
    coeff = charge_C / (4.0 * math.pi * EPSILON_0)

    if radius == 0.0:
        phi = coeff * math.sqrt(2.0 / math.pi) / sigma_m
        return phi, np.zeros(3)

    u = radius / (math.sqrt(2.0) * sigma_m)
    exp_term = math.exp(-0.5 * radius2 / sigma2)
    erf_u = math.erf(u)

    phi = coeff * erf_u / radius
    shape = erf_u - math.sqrt(2.0 / math.pi) * (radius / sigma_m) * exp_term
    e_prime = coeff * shape / (radius**3) * r_prime_m
    return phi, e_prime


def source_rest_coordinate(position_m: np.ndarray, time_s: float, source: Source) -> np.ndarray:
    """Map a lab event position at fixed lab time into the source rest frame."""
    displacement = position_m - source.center(time_s)
    parallel, perpendicular = split_parallel_perpendicular(displacement, source.beta_hat)
    return perpendicular + source.gamma * parallel


def boost_rest_fields_to_lab(e_prime: np.ndarray, source: Source) -> tuple[np.ndarray, np.ndarray]:
    """Boost source-rest fields with B'=0 into the lab frame."""
    parallel, perpendicular = split_parallel_perpendicular(e_prime, source.beta_hat)
    e_lab = parallel + source.gamma * perpendicular
    b_lab = np.cross(source.beta * C_LIGHT, e_lab) / (C_LIGHT * C_LIGHT)
    return e_lab, b_lab


def source_lab_fields(position_m: np.ndarray, time_s: float, source: Source) -> tuple[float, np.ndarray, np.ndarray]:
    """Return phi', E_lab, B_lab from one uniformly moving Gaussian source."""
    r_prime = source_rest_coordinate(position_m, time_s, source)
    phi_prime, e_prime = gaussian_rest_field(r_prime, source.charge_C, source.sigma_m)
    e_lab, b_lab = boost_rest_fields_to_lab(e_prime, source)
    return phi_prime, e_lab, b_lab


def total_lab_fields(position_m: np.ndarray, time_s: float, sources: Iterable[Source]) -> tuple[np.ndarray, np.ndarray]:
    """Return total lab-frame E and B at one witness event."""
    e_total = np.zeros(3)
    b_total = np.zeros(3)
    for source in sources:
        _phi_prime, e_lab, b_lab = source_lab_fields(position_m, time_s, source)
        e_total += e_lab
        b_total += b_lab
    return e_total, b_total


def boris_kick_opalx_units(
    momentum: np.ndarray,
    e_field: np.ndarray,
    b_field: np.ndarray,
    dt_s: float,
    charge_units: float = ELECTRON_CHARGE_UNITS,
    mass_eV: float = ELECTRON_REST_EV,
) -> np.ndarray:
    """Apply the OPALX Boris kick to normalized momentum P=beta*gamma."""
    p = momentum.astype(float).copy()
    p += 0.5 * dt_s * charge_units * C_LIGHT / mass_eV * e_field

    gamma = math.sqrt(1.0 + float(np.dot(p, p)))
    t_vec = 0.5 * dt_s * charge_units * C_LIGHT * C_LIGHT / (gamma * mass_eV) * b_field
    w_vec = p + np.cross(p, t_vec)
    s_vec = 2.0 / (1.0 + float(np.dot(t_vec, t_vec))) * t_vec
    p += np.cross(w_vec, s_vec)

    p += 0.5 * dt_s * charge_units * C_LIGHT / mass_eV * e_field
    return p


def advance_witness_boris(
    position_m: np.ndarray,
    momentum: np.ndarray,
    charge_units: float,
    time_s: float,
    dt_s: float,
    sources: Iterable[Source],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Advance one witness by one step and return the fields used."""
    e_field, b_field = total_lab_fields(position_m, time_s, sources)
    p_new = boris_kick_opalx_units(momentum, e_field, b_field, dt_s, charge_units)

    # Trapezoidal drift keeps the manufactured trajectory symmetric under
    # timestep refinement while preserving OPALX's normalized momentum units.
    v_old = velocity_from_momentum(momentum)
    v_new = velocity_from_momentum(p_new)
    r_new = position_m + 0.5 * (v_old + v_new) * dt_s
    return r_new, p_new, e_field, b_field


def integrate_witness_particle(
    particle: WitnessParticle,
    sources: tuple[Source, ...],
    duration_s: float,
    dt_s: float,
    t0_s: float = 0.0,
) -> pd.DataFrame:
    """Integrate one witness particle and return a trajectory dataframe."""
    if duration_s <= 0.0:
        raise ValueError("duration_s must be positive")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")

    steps = int(math.ceil(duration_s / dt_s))
    position = particle.position_m.astype(float).copy()
    momentum = particle.momentum.astype(float).copy()
    rows = []

    for step in range(steps + 1):
        time_s = t0_s + min(step * dt_s, duration_s)
        e_field, b_field = total_lab_fields(position, time_s, sources)
        gamma = gamma_from_momentum(momentum)
        rows.append(
            {
                "species": particle.name,
                "step": step,
                "t_ps": 1.0e12 * (time_s - t0_s),
                "x_m": position[0],
                "y_m": position[1],
                "z_m": position[2],
                "px": momentum[0],
                "py": momentum[1],
                "pz": momentum[2],
                "gamma": gamma,
                "kinetic_keV": (gamma - 1.0) * ELECTRON_REST_EV / 1.0e3,
                "Ex_V_per_m": e_field[0],
                "Ey_V_per_m": e_field[1],
                "Ez_V_per_m": e_field[2],
                "Bx_T": b_field[0],
                "By_T": b_field[1],
                "Bz_T": b_field[2],
                "E_abs_V_per_m": float(np.linalg.norm(e_field)),
                "B_abs_T": float(np.linalg.norm(b_field)),
            }
        )
        if step == steps:
            break
        this_dt = min(dt_s, t0_s + duration_s - time_s)
        position, momentum, _e_field, _b_field = advance_witness_boris(
            position, momentum, particle.charge_units, time_s, this_dt, sources
        )

    return pd.DataFrame(rows)


def trajectory_row(
    particle_name: str,
    step: int,
    time_s: float,
    t0_s: float,
    position: np.ndarray,
    momentum: np.ndarray,
    sources: Iterable[Source],
) -> dict[str, float | int | str]:
    """Return one sampled trajectory row."""
    e_field, b_field = total_lab_fields(position, time_s, sources)
    gamma = gamma_from_momentum(momentum)
    return {
        "species": particle_name,
        "step": step,
        "t_ps": 1.0e12 * (time_s - t0_s),
        "x_m": position[0],
        "y_m": position[1],
        "z_m": position[2],
        "px": momentum[0],
        "py": momentum[1],
        "pz": momentum[2],
        "gamma": gamma,
        "kinetic_keV": (gamma - 1.0) * ELECTRON_REST_EV / 1.0e3,
        "Ex_V_per_m": e_field[0],
        "Ey_V_per_m": e_field[1],
        "Ez_V_per_m": e_field[2],
        "Bx_T": b_field[0],
        "By_T": b_field[1],
        "Bz_T": b_field[2],
        "E_abs_V_per_m": float(np.linalg.norm(e_field)),
        "B_abs_T": float(np.linalg.norm(b_field)),
    }


def integrate_witness_particle_multirate(
    particle: WitnessParticle,
    sources: tuple[Source, ...],
    duration_s: float,
    fine_duration_s: float,
    fine_dt_s: float,
    tail_dt_s: float,
    output_dt_s: float,
    t0_s: float = 0.0,
) -> pd.DataFrame:
    """Integrate one witness with fine overlap steps and sparse output rows."""
    if duration_s <= 0.0:
        raise ValueError("duration_s must be positive")
    if fine_duration_s < 0.0:
        raise ValueError("fine_duration_s must be non-negative")
    if fine_dt_s <= 0.0 or tail_dt_s <= 0.0 or output_dt_s <= 0.0:
        raise ValueError("all timesteps must be positive")

    position = particle.position_m.astype(float).copy()
    momentum = particle.momentum.astype(float).copy()
    rows = []
    step = 0
    time_s = t0_s
    next_output_s = t0_s
    fine_end_s = t0_s + min(fine_duration_s, duration_s)
    end_s = t0_s + duration_s
    tolerance_s = 1.0e-18

    rows.append(trajectory_row(particle.name, step, time_s, t0_s, position, momentum, sources))
    next_output_s += output_dt_s

    while time_s < end_s - tolerance_s:
        dt_s = fine_dt_s if time_s < fine_end_s - tolerance_s else tail_dt_s
        boundary_s = fine_end_s if time_s < fine_end_s - tolerance_s else end_s
        this_dt = min(dt_s, boundary_s - time_s, end_s - time_s)
        position, momentum, _e_field, _b_field = advance_witness_boris(
            position, momentum, particle.charge_units, time_s, this_dt, sources
        )
        time_s += this_dt
        step += 1
        if time_s + tolerance_s >= next_output_s or time_s >= end_s - tolerance_s:
            rows.append(trajectory_row(particle.name, step, time_s, t0_s, position, momentum, sources))
            while next_output_s <= time_s + tolerance_s:
                next_output_s += output_dt_s

    return pd.DataFrame(rows)


def make_default_sources(args: argparse.Namespace) -> tuple[Source, Source, Source]:
    """Build symmetric and weakened electron-bunch sources for the table."""
    beta = beta_from_kinetic_energy(args.source_kinetic_MeV)
    z_sep = args.half_separation_m
    base = {
        "sigma_m": args.sigma_m,
        "t0_s": 0.0,
    }
    left = Source(
        "left",
        args.charge_C,
        beta=np.array([0.0, 0.0, +beta]),
        center_t0_m=np.array([0.0, 0.0, -z_sep]),
        **base,
    )
    right = Source(
        "right",
        args.charge_C,
        beta=np.array([0.0, 0.0, -beta]),
        center_t0_m=np.array([0.0, 0.0, +z_sep]),
        **base,
    )
    weak_right = Source(
        "right_80pct",
        0.8 * args.charge_C,
        beta=np.array([0.0, 0.0, -beta]),
        center_t0_m=np.array([0.0, 0.0, +z_sep]),
        **base,
    )
    return left, right, weak_right


def default_cases(args: argparse.Namespace) -> list[WitnessCase]:
    """Return the manufactured validation cases used in the LaTeX note."""
    left, right, weak_right = make_default_sources(args)
    p0 = np.array(args.witness_momentum, dtype=float)
    x0 = args.witness_x_m
    return [
        WitnessCase(
            "symmetric center",
            "identical head-on sources, witness at IP",
            np.array([0.0, 0.0, 0.0]),
            p0,
            (left, right),
        ),
        WitnessCase(
            "symmetric off-axis",
            "identical head-on sources, transverse witness offset",
            np.array([x0, 0.0, 0.0]),
            p0,
            (left, right),
        ),
        WitnessCase(
            "single boosted source",
            "one source only, off-axis witness in the source-center plane",
            np.array([x0, 0.0, 0.0]),
            p0,
            (left,),
        ),
        WitnessCase(
            "asymmetric collision",
            "right source charge reduced to 80 percent",
            np.array([x0, 0.0, 0.0]),
            p0,
            (left, weak_right),
        ),
    ]


def compute_case(case: WitnessCase, time_s: float, dt_s: float) -> dict[str, object]:
    """Evaluate one manufactured case."""
    e_field, b_field = total_lab_fields(case.position_m, time_s, case.sources)
    p_after = boris_kick_opalx_units(case.momentum, e_field, b_field, dt_s)
    dp = p_after - case.momentum
    return {
        "case": case,
        "E": e_field,
        "B": b_field,
        "P_after": p_after,
        "dP": dp,
        "E_norm": float(np.linalg.norm(e_field)),
        "B_norm": float(np.linalg.norm(b_field)),
        "dP_norm": float(np.linalg.norm(dp)),
    }


DEFAULT_LATEX_TABLE = NOTE_DIR / "boosted_gaussian_witness_initial_cases_table.tex"
DEFAULT_PAIR_KINEMATICS_TABLE = NOTE_DIR / "boosted_gaussian_witness_pair_kinematics_table.tex"


def fmt_vec(vector: np.ndarray, scale: float = 1.0) -> str:
    """Format a three-vector compactly."""
    return "(" + ", ".join(f"{scale * value:.3e}" for value in vector) + ")"


def fmt_table_scalar(value: float, precision: int = 3) -> str:
    """Format table values in fixed notation when the chosen units are readable."""
    if not math.isfinite(value):
        return str(value)
    if abs(value) < 0.5 * 10.0 ** (-precision):
        value = 0.0
    if abs(value) < 1.0e4:
        return f"{value:.{precision}f}"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / 10.0**exponent
    return rf"{mantissa:.{precision}f}\times 10^{{{exponent:d}}}"


def fmt_table_vec(vector: np.ndarray, scale: float = 1.0, precision: int = 3) -> str:
    """Format a three-vector for generated LaTeX tables."""
    return "(" + ", ".join(fmt_table_scalar(scale * value, precision) for value in vector) + ")"


def print_text_table(rows: list[dict[str, object]]) -> None:
    """Print a compact text table."""
    header = (
        f"{'case':24s} {'r_w [mm]':24s} {'|E| [GV/m]':>12s} "
        f"{'|B| [T]':>10s} {'Delta P':>26s}"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        case = row["case"]
        assert isinstance(case, WitnessCase)
        print(
            f"{case.name:24s} {fmt_table_vec(case.position_m, 1.0e3):24s} "
            f"{fmt_table_scalar(1.0e-9 * row['E_norm']):>12s} "
            f"{fmt_table_scalar(row['B_norm']):>10s} "
            f"{fmt_table_vec(row['dP']):>26s}"
        )


def latex_table(rows: list[dict[str, object]]) -> str:
    """Return a LaTeX tabular containing the manufactured table."""
    lines = [
        r"% Generated by boosted_gaussian_witness.py. Do not edit by hand.",
        r"\begin{tabular}{lllll}",
        r"\hline",
        r"case & $\mathbf r_w$ [mm] & $|\mathbf E|$ [GV/m] & $|\mathbf B|$ [T] & $\Delta\mathbf P$ \\",
        r"\hline",
    ]
    for row in rows:
        case = row["case"]
        assert isinstance(case, WitnessCase)
        lines.append(
            f"{case.name} & "
            f"${fmt_table_vec(case.position_m, 1.0e3)}$ & "
            f"${fmt_table_scalar(1.0e-9 * row['E_norm'])}$ & "
            f"${fmt_table_scalar(row['B_norm'])}$ & "
            f"${fmt_table_vec(row['dP'])}$ \\\\"
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return "\n".join(lines)


def latex_pair_kinematics_table(
    trajectories: pd.DataFrame,
    field_scenarios: tuple[tuple[str, str, tuple[Source, ...]], ...],
) -> str:
    """Return a LaTeX table summarizing field-on pair kinematics against free drift."""
    free_final = trajectories[trajectories["scenario"] == "no_fields"].set_index("species")
    lines = [
        r"% Generated by boosted_gaussian_witness.py. Do not edit by hand.",
        r"\begin{tabular}{lllllll}",
        r"\hline",
        (
            r"case & particle & $\mathbf r_f$ [mm] & $\mathbf P_f$ "
            r"& $K_f$ [keV] & $\Delta\mathbf r$ [$\mu$m] & $\Delta K$ [eV] \\"
        ),
        r"\hline",
    ]

    for scenario, label, _sources in field_scenarios:
        field_final = trajectories[trajectories["scenario"] == scenario].set_index("species")
        for species, particle_label in (("electron", r"e^-"), ("positron", r"e^+")):
            field_row = field_final.loc[species].iloc[-1]
            free_row = free_final.loc[species].iloc[-1]
            final_r = field_row[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
            final_p = field_row[["px", "py", "pz"]].to_numpy(dtype=float)
            delta_r = (
                field_row[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
                - free_row[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
            )
            delta_k_eV = 1.0e3 * (field_row["kinetic_keV"] - free_row["kinetic_keV"])
            lines.append(
                f"{label} & "
                f"${particle_label}$ & "
                f"${fmt_table_vec(final_r, 1.0e3, 6)}$ & "
                f"${fmt_table_vec(final_p, 1.0, 6)}$ & "
                f"${fmt_table_scalar(field_row['kinetic_keV'], 3)}$ & "
                f"${fmt_table_vec(delta_r, 1.0e6, 3)}$ & "
                f"${fmt_table_scalar(delta_k_eV, 3)}$ \\\\"
            )

    lines.extend([r"\hline", r"\end{tabular}"])
    return "\n".join(lines)


def configure_plot_environment(output_prefix: Path) -> None:
    """Use writable local cache directories for headless Matplotlib rendering."""
    cache_dir = output_prefix.parent / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))


def make_breit_wheeler_pair(args: argparse.Namespace) -> tuple[WitnessParticle, WitnessParticle]:
    """Return an electron/positron pair born back-to-back at the IP."""
    direction = np.asarray(args.pair_direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        raise ValueError("--pair-direction must not be the zero vector")
    direction /= norm

    p_abs = momentum_magnitude_from_kinetic_energy_keV(args.pair_energy_keV)
    origin = np.asarray(args.pair_origin_m, dtype=float)
    electron = WitnessParticle(
        name="electron",
        charge_units=-1.0,
        position_m=origin.copy(),
        momentum=p_abs * direction,
    )
    positron = WitnessParticle(
        name="positron",
        charge_units=+1.0,
        position_m=origin.copy(),
        momentum=-p_abs * direction,
    )
    return electron, positron


def field_magnitude_plane(
    sources: tuple[Source, ...],
    time_s: float,
    plane: str,
    span_m: float,
    n_grid: int,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Evaluate |E| on an xy or xz plane through the IP."""
    axis = np.linspace(-0.5 * span_m, 0.5 * span_m, n_grid)
    values = np.zeros((n_grid, n_grid), dtype=float)
    for i, a in enumerate(axis):
        for j, b in enumerate(axis):
            if plane == "xy":
                position = np.array([b, a, 0.0])
            elif plane == "xz":
                position = np.array([b, 0.0, a])
            else:
                raise ValueError(f"unsupported plane {plane!r}")
            e_field, _b_field = total_lab_fields(position, time_s, sources)
            values[i, j] = float(np.linalg.norm(e_field))
    extent_mm = (
        -0.5 * span_m * 1.0e3,
        0.5 * span_m * 1.0e3,
        -0.5 * span_m * 1.0e3,
        0.5 * span_m * 1.0e3,
    )
    return values, extent_mm


def field_magnitude_volume(
    sources: tuple[Source, ...],
    time_s: float,
    span_m: float,
    n_grid: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate |E| on a cubic grid for a compact 3D visualization."""
    axis = np.linspace(-0.5 * span_m, 0.5 * span_m, n_grid)
    points = []
    values = []
    for z in axis:
        for y in axis:
            for x in axis:
                position = np.array([x, y, z])
                e_field, _b_field = total_lab_fields(position, time_s, sources)
                points.append((1.0e3 * x, 1.0e3 * y, 1.0e3 * z))
                values.append(float(np.linalg.norm(e_field)))
    xyz = np.asarray(points, dtype=float)
    return xyz[:, 0], xyz[:, 1], xyz[:, 2], np.asarray(values, dtype=float)


def log_color_limits(*arrays: np.ndarray) -> tuple[float, float]:
    """Return robust log10 color limits for positive field magnitudes."""
    positive = np.concatenate([array[np.isfinite(array) & (array > 0.0)] for array in arrays])
    if positive.size == 0:
        return -1.0, 1.0
    vmin = float(np.percentile(positive, 1.0))
    vmax = float(np.percentile(positive, 99.5))
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin >= vmax:
        vmin = float(np.min(positive))
        vmax = float(np.max(positive))
    if vmin <= 0.0 or vmin >= vmax:
        vmax = max(vmax, 1.0)
        vmin = 1.0e-12 * vmax
    return math.log10(vmin), math.log10(vmax)


def log_field(values: np.ndarray) -> np.ndarray:
    """Return log10(|E|), with exact nulls masked for plotting."""
    logged = np.full_like(values, np.nan, dtype=float)
    mask = values > 0.0
    logged[mask] = np.log10(values[mask])
    return logged


def plot_field_magnitude_t0(
    field_scenarios: tuple[tuple[str, str, tuple[Source, ...]], ...],
    output: Path,
    span_m: float,
    n_grid: int,
) -> None:
    """Plot |E| at t=0 in two central slices and a sampled 3D volume."""
    configure_plot_environment(output)
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 220, "font.size": 11})
    volume_grid = max(17, min(35, n_grid))
    scenario_data = []
    color_arrays = []
    for scenario_name, scenario_label, sources in field_scenarios:
        xy_values, xy_extent = field_magnitude_plane(sources, 0.0, "xy", span_m, n_grid)
        xz_values, xz_extent = field_magnitude_plane(sources, 0.0, "xz", span_m, n_grid)
        x_mm, y_mm, z_mm, volume_values = field_magnitude_volume(sources, 0.0, span_m, volume_grid)
        scenario_data.append(
            {
                "name": scenario_name,
                "label": scenario_label,
                "xy_values": xy_values,
                "xy_extent": xy_extent,
                "xz_values": xz_values,
                "xz_extent": xz_extent,
                "x_mm": x_mm,
                "y_mm": y_mm,
                "z_mm": z_mm,
                "volume_values": volume_values,
            }
        )
        color_arrays.extend([xy_values, xz_values, volume_values])

    log_vmin, log_vmax = log_color_limits(*color_arrays)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#111111")

    n_rows = len(scenario_data)
    fig = plt.figure(figsize=(16, 5.2 * n_rows), constrained_layout=True)
    axes_for_colorbar = []
    mappable = None
    for row, data in enumerate(scenario_data):
        ax_xy = fig.add_subplot(n_rows, 3, 3 * row + 1)
        ax_xz = fig.add_subplot(n_rows, 3, 3 * row + 2)
        ax_3d = fig.add_subplot(n_rows, 3, 3 * row + 3, projection="3d")
        axes_for_colorbar.extend([ax_xy, ax_xz, ax_3d])
        panels = (
            (
                ax_xy,
                data["xy_values"],
                data["xy_extent"],
                rf"{data['label']}: $|\mathbf{{E}}|$ at $z=0$",
                "x [mm]",
                "y [mm]",
            ),
            (
                ax_xz,
                data["xz_values"],
                data["xz_extent"],
                rf"{data['label']}: $|\mathbf{{E}}|$ at $y=0$",
                "x [mm]",
                "z [mm]",
            ),
        )
        for ax, values, extent, title, xlabel, ylabel in panels:
            image = ax.imshow(
                log_field(values),
                origin="lower",
                extent=extent,
                aspect="equal",
                cmap=cmap,
                vmin=log_vmin,
                vmax=log_vmax,
            )
            mappable = image
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.plot(0.0, 0.0, marker="+", markersize=9, color="white", mew=1.5)

        volume_values = data["volume_values"]
        finite_volume = volume_values[np.isfinite(volume_values) & (volume_values > 0.0)]
        if finite_volume.size:
            threshold = float(np.percentile(finite_volume, 75.0))
            mask = volume_values >= threshold
            if np.count_nonzero(mask) > 7000:
                threshold = float(np.partition(volume_values, -7000)[-7000])
                mask = volume_values >= threshold
            mappable = ax_3d.scatter(
                data["x_mm"][mask],
                data["y_mm"][mask],
                data["z_mm"][mask],
                c=log_field(volume_values[mask]),
                cmap=cmap,
                vmin=log_vmin,
                vmax=log_vmax,
                s=5,
                alpha=0.55,
                linewidths=0,
            )
        ax_3d.scatter([0.0], [0.0], [0.0], marker="+", s=60, color="black", linewidths=1.5)
        ax_3d.set_title(rf"{data['label']}: 3D sampled $|\mathbf{{E}}|$")
        ax_3d.set_xlabel("x [mm]")
        ax_3d.set_ylabel("y [mm]")
        ax_3d.set_zlabel("z [mm]")
        set_axes_equal(ax_3d)
        ax_3d.view_init(elev=22, azim=-45)

    assert mappable is not None
    colorbar = fig.colorbar(mappable, ax=axes_for_colorbar, shrink=0.82)
    colorbar.set_label(r"$\log_{10}(|\mathbf{E}| / [\mathrm{V/m}])$")
    fig.suptitle("Boosted Gaussian source field magnitude at pair creation")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def set_axes_equal(ax) -> None:
    """Set equal scale on a Matplotlib 3D axis."""
    limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()], dtype=float)
    centers = limits.mean(axis=1)
    radius = 0.5 * np.max(limits[:, 1] - limits[:, 0])
    for center, setter in zip(centers, (ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d)):
        setter(center - radius, center + radius)


def expand_flat_axis(ax, axis: str, values: np.ndarray, reference_span: float) -> None:
    """Expand a nearly flat 2D axis so planar trajectories remain legible."""
    center = 0.5 * (float(np.min(values)) + float(np.max(values)))
    span = float(np.max(values) - np.min(values))
    if span < 1.0e-9:
        span = max(0.05 * reference_span, 0.05)
    else:
        span *= 1.15
    limits = (center - 0.5 * span, center + 0.5 * span)
    if axis == "x":
        ax.set_xlim(*limits)
    elif axis == "y":
        ax.set_ylim(*limits)
    else:
        raise ValueError(f"unsupported axis {axis!r}")


def plot_pair_paths(
    trajectories: pd.DataFrame,
    sources: tuple[Source, ...],
    output: Path,
    beambeam_radius_m: float,
    beambeam_length_m: float,
) -> None:
    """Plot pair trajectories and basic energy diagnostics."""
    configure_plot_environment(output)
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 220, "font.size": 11})
    colors = {"electron": "#2166ac", "positron": "#b2182b"}
    linestyles = {"symmetric_collision": "-", "asymmetric_collision": "-.", "no_fields": "--"}
    alphas = {"symmetric_collision": 1.0, "asymmetric_collision": 0.9, "no_fields": 0.55}
    scenario_labels = {
        "symmetric_collision": "symmetric collision",
        "asymmetric_collision": "asymmetric collision",
        "no_fields": "no fields",
    }
    fig = plt.figure(figsize=(13.0, 8.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax3d = fig.add_subplot(grid[0, 0], projection="3d")
    axxz = fig.add_subplot(grid[0, 1])
    axk = fig.add_subplot(grid[1, 0])
    axdelta = fig.add_subplot(grid[1, 1])
    all_x_mm = []
    all_z_um = []
    species_handles = {}
    scenario_handles = {}

    group_columns = ["scenario", "species"] if "scenario" in trajectories.columns else ["species"]
    for group_key, frame in trajectories.groupby(group_columns, sort=False):
        if len(group_columns) == 2:
            scenario, species = group_key
        else:
            species = group_key
            scenario = "symmetric_collision"
        color = colors.get(species, None)
        linestyle = linestyles.get(scenario, "-")
        alpha = alphas.get(scenario, 1.0)
        label = r"$e^-$" if species == "electron" else r"$e^+$"
        x_mm = 1.0e3 * frame["x_m"].to_numpy()
        y_mm = 1.0e3 * frame["y_m"].to_numpy()
        z_mm = 1.0e3 * frame["z_m"].to_numpy()
        z_um = 1.0e6 * frame["z_m"].to_numpy()
        all_x_mm.append(x_mm)
        all_z_um.append(z_um)
        line3d, = ax3d.plot(x_mm, y_mm, z_mm, color=color, ls=linestyle, lw=2.2, alpha=alpha)
        ax3d.scatter(x_mm[0], y_mm[0], z_mm[0], color=color, marker="o", s=28, alpha=alpha)
        ax3d.scatter(x_mm[-1], y_mm[-1], z_mm[-1], color=color, marker="s", s=28, alpha=alpha)
        linexz, = axxz.plot(x_mm, z_um, color=color, ls=linestyle, lw=2.2, alpha=alpha)
        axxz.scatter(x_mm[0], z_um[0], color=color, marker="o", s=28, alpha=alpha)
        axxz.scatter(x_mm[-1], z_um[-1], color=color, marker="s", s=28, alpha=alpha)
        axk.plot(
            frame["t_ps"].to_numpy(),
            1.0e3 * (frame["kinetic_keV"].to_numpy() - frame["kinetic_keV"].to_numpy()[0]),
            color=color,
            ls=linestyle,
            lw=2.0,
            alpha=alpha,
        )
        species_handles.setdefault(species, line3d)
        scenario_handles.setdefault(scenario, linexz)

    time_ps = trajectories["t_ps"].to_numpy()
    source_times = np.linspace(time_ps.min(), time_ps.max(), 200) * 1.0e-12
    for source in sources:
        centers = np.asarray([source.center(t) for t in source_times])
        ax3d.plot(
            1.0e3 * centers[:, 0],
            1.0e3 * centers[:, 1],
            1.0e3 * centers[:, 2],
            color="#7f7f7f",
            ls="--",
            lw=1.0,
            alpha=0.7,
        )
        axxz.plot(
            1.0e3 * centers[:, 0],
            1.0e6 * centers[:, 2],
            color="#7f7f7f",
            ls="--",
            lw=1.0,
            alpha=0.6,
        )

    ax3d.set_title("3D witness trajectories")
    ax3d.set_xlabel("x [mm]")
    ax3d.set_ylabel("y [mm]")
    ax3d.set_zlabel("z [mm]")
    species_names = [name for name in ("electron", "positron") if name in species_handles]
    species_legend = ax3d.legend(
        [species_handles[name] for name in species_names],
        [r"$e^-$" if name == "electron" else r"$e^+$" for name in species_names],
        frameon=False,
        loc="upper right",
    )
    ax3d.add_artist(species_legend)
    set_axes_equal(ax3d)
    ax3d.view_init(elev=22, azim=-58)

    axxz.set_title("x-z projection")
    axxz.set_xlabel("x [mm]")
    axxz.set_ylabel(r"z [$\mu$m]")
    axxz.grid(True, alpha=0.25)
    radius_mm = 1.0e3 * beambeam_radius_m
    half_length_um = 0.5e6 * beambeam_length_m
    axxz.axvline(radius_mm, color="0.35", ls=":", lw=1.1, label="15 cm radius")
    axxz.axvline(-radius_mm, color="0.35", ls=":", lw=1.1)
    axxz.axhspan(-half_length_um, half_length_um, color="0.92", zorder=-4, label="32 cm window")
    all_x = np.concatenate(all_x_mm)
    all_z = np.concatenate(all_z_um)
    x_span = max(float(np.max(all_x) - np.min(all_x)), 1.0)
    expand_flat_axis(axxz, "x", all_x, x_span)
    expand_flat_axis(axxz, "y", all_z, x_span)
    axxz.legend(
        [
            scenario_handles[name]
            for name in ("symmetric_collision", "asymmetric_collision", "no_fields")
            if name in scenario_handles
        ],
        [
            scenario_labels[name]
            for name in ("symmetric_collision", "asymmetric_collision", "no_fields")
            if name in scenario_handles
        ],
        frameon=False,
        fontsize=10,
        loc="upper right",
    )

    axk.axhline(0.0, color="#666666", lw=0.8, alpha=0.6)
    axk.set_title("kinetic energy change")
    axk.set_xlabel("t [ps]")
    axk.set_ylabel(r"$\Delta K$ [eV]")
    axk.grid(True, alpha=0.25)

    if "scenario" in trajectories.columns:
        for scenario in ("symmetric_collision", "asymmetric_collision"):
            for species in ("electron", "positron"):
                field_frame = trajectories[
                    (trajectories["scenario"] == scenario) & (trajectories["species"] == species)
                ].sort_values("step")
                free_frame = trajectories[
                    (trajectories["scenario"] == "no_fields") & (trajectories["species"] == species)
                ].sort_values("step")
                if field_frame.empty or free_frame.empty:
                    continue
                delta_um = 1.0e6 * (
                    field_frame[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
                    - free_frame[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
                )
                species_label = r"$e^-$" if species == "electron" else r"$e^+$"
                axdelta.plot(
                    field_frame["t_ps"].to_numpy(),
                    delta_um[:, 0],
                    color=colors.get(species),
                    ls=linestyles.get(scenario),
                    alpha=alphas.get(scenario, 1.0),
                    lw=2.0,
                    label=f"{species_label}, {scenario_labels[scenario]}",
                )
    axdelta.axhline(0.0, color="#666666", lw=0.8, alpha=0.6)
    axdelta.set_title("field-on minus free propagation")
    axdelta.set_xlabel("t [ps]")
    axdelta.set_ylabel(r"$\Delta x$ [$\mu$m]")
    axdelta.grid(True, alpha=0.25)
    axdelta.legend(frameon=False, fontsize=9)
    fig.suptitle(
            "Breit-Wheeler e-/e+ pair in boosted Gaussian source fields "
            f"(BeamBeam radius {beambeam_radius_m * 1.0e2:.0f} cm, "
            f"length {beambeam_length_m * 1.0e2:.0f} cm)")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def run_pair_demo(args: argparse.Namespace) -> int:
    """Integrate a back-to-back e-/e+ pair and write plots/data."""
    left, right, weak_right = make_default_sources(args)
    symmetric_sources = (left, right)
    asymmetric_sources = (left, weak_right)
    electron, positron = make_breit_wheeler_pair(args)
    duration_s = args.pair_duration_ps * 1.0e-12
    dt_s = args.pair_dt_ps * 1.0e-12
    fine_dt_s = None if args.pair_fine_dt_ps is None else args.pair_fine_dt_ps * 1.0e-12
    fine_duration_s = args.pair_fine_duration_ps * 1.0e-12
    output_dt_s = args.pair_output_dt_ps * 1.0e-12

    trajectory_frames = []
    scenarios = (
        ("symmetric_collision", symmetric_sources),
        ("asymmetric_collision", asymmetric_sources),
        ("no_fields", ()),
    )
    for scenario, scenario_sources in scenarios:
        for particle in (electron, positron):
            if fine_dt_s is None:
                frame = integrate_witness_particle(particle, scenario_sources, duration_s, dt_s)
            else:
                frame = integrate_witness_particle_multirate(
                    particle,
                    scenario_sources,
                    duration_s,
                    fine_duration_s,
                    fine_dt_s,
                    dt_s,
                    output_dt_s,
                )
            frame.insert(0, "scenario", scenario)
            trajectory_frames.append(frame)
    trajectories = pd.concat(trajectory_frames, ignore_index=True)

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_prefix.with_name(f"{output_prefix.name}_trajectories.csv")
    field_path = output_prefix.with_name(f"{output_prefix.name}_field_t0.png")
    path_path = output_prefix.with_name(f"{output_prefix.name}_paths.png")
    trajectories.to_csv(csv_path, index=False)
    field_scenarios = (
        ("symmetric_collision", "symmetric collision", symmetric_sources),
        ("asymmetric_collision", "asymmetric collision", asymmetric_sources),
    )
    pair_table_path = args.write_pair_latex_table
    pair_table_path.parent.mkdir(parents=True, exist_ok=True)
    pair_table_path.write_text(
        latex_pair_kinematics_table(trajectories, field_scenarios) + "\n",
        encoding="utf-8",
    )
    plot_field_magnitude_t0(field_scenarios, field_path, args.field_span_mm * 1.0e-3, args.field_grid)
    plot_pair_paths(trajectories, symmetric_sources, path_path, args.beambeam_radius_m, args.beambeam_length_m)

    p_abs = momentum_magnitude_from_kinetic_energy_keV(args.pair_energy_keV)
    beta_pair = p_abs / math.sqrt(1.0 + p_abs * p_abs)
    print("Breit-Wheeler pair demo")
    print(f"pair kinetic energy = {args.pair_energy_keV:.6g} keV")
    print(f"|P| = beta*gamma = {p_abs:.12g}, beta = {beta_pair:.12g}")
    if fine_dt_s is None:
        print(f"duration = {args.pair_duration_ps:.6g} ps, dt = {args.pair_dt_ps:.6g} ps")
    else:
        print(
            f"duration = {args.pair_duration_ps:.6g} ps, "
            f"fine dt = {args.pair_fine_dt_ps:.6g} ps to {args.pair_fine_duration_ps:.6g} ps, "
            f"tail dt = {args.pair_dt_ps:.6g} ps, output dt = {args.pair_output_dt_ps:.6g} ps"
        )
    for scenario, _label, scenario_sources in field_scenarios:
        e0, b0 = total_lab_fields(np.asarray(args.pair_origin_m, dtype=float), 0.0, scenario_sources)
        print(f"{scenario:20s} E(t=0, origin) = {fmt_vec(e0)} V/m, |E|={np.linalg.norm(e0):.6e} V/m")
        print(f"{scenario:20s} B(t=0, origin) = {fmt_vec(b0)} T, |B|={np.linalg.norm(b0):.6e} T")
    for (scenario, species), frame in trajectories.groupby(["scenario", "species"], sort=False):
        final = frame.iloc[-1]
        final_r = final[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
        final_p = final[["px", "py", "pz"]].to_numpy(dtype=float)
        print(
            f"{scenario:16s} {species:8s} final r [mm]={fmt_vec(final_r, 1.0e3)} "
            f"P={fmt_vec(final_p)} K={final['kinetic_keV']:.6f} keV"
        )
        print(
            f"{scenario:16s} {species:8s} max |E|={frame['E_abs_V_per_m'].max():.6e} V/m, "
            f"max |B|={frame['B_abs_T'].max():.6e} T"
        )
    free_final = trajectories[trajectories["scenario"] == "no_fields"].set_index("species")
    for scenario, _label, _sources in field_scenarios:
        field_final = trajectories[trajectories["scenario"] == scenario].set_index("species")
        for species in ("electron", "positron"):
            field_row = field_final.loc[species].iloc[-1]
            free_row = free_final.loc[species].iloc[-1]
            dr_mm = 1.0e3 * (
                field_row[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
                - free_row[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)
            )
            dk_eV = 1.0e3 * (field_row["kinetic_keV"] - free_row["kinetic_keV"])
            print(f"{scenario:20s} minus free {species:8s}: dr [mm]={fmt_vec(dr_mm)}, dK={dk_eV:.6e} eV")
    print(f"wrote {csv_path}")
    print(f"wrote {pair_table_path}")
    print(f"wrote {field_path}")
    print(f"wrote {path_path}")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Lorentz-transformed Gaussian manufactured witness-kick cases."
    )
    parser.add_argument("--sigma-m", type=float, default=0.6e-3, help="rest-frame Gaussian sigma")
    parser.add_argument(
        "--charge-C",
        type=float,
        default=DEFAULT_SOURCE_CHARGE_C,
        help="source bunch charge",
    )
    parser.add_argument(
        "--source-kinetic-MeV",
        type=float,
        default=245.0,
        help="source electron kinetic energy used to set beta",
    )
    parser.add_argument(
        "--half-separation-m",
        type=float,
        default=0.0,
        help="lab-frame source half separation at t=0",
    )
    parser.add_argument("--witness-x-m", type=float, default=5.0e-4)
    parser.add_argument("--time-s", type=float, default=0.0)
    parser.add_argument("--dt-s", type=float, default=1.0e-12)
    parser.add_argument(
        "--witness-momentum",
        type=float,
        nargs=3,
        default=(0.0, 0.0, 0.0),
        metavar=("PX", "PY", "PZ"),
        help="initial witness normalized momentum beta*gamma",
    )
    parser.add_argument("--latex-table", action="store_true", help="print a LaTeX table")
    parser.add_argument(
        "--write-latex-table",
        type=Path,
        default=DEFAULT_LATEX_TABLE,
        help="write the LaTeX table to a file",
    )
    parser.add_argument(
        "--pair-demo",
        action="store_true",
        help="integrate a back-to-back 313 keV e-/e+ pair and write plots",
    )
    parser.add_argument(
        "--write-pair-latex-table",
        type=Path,
        default=DEFAULT_PAIR_KINEMATICS_TABLE,
        help="write the pair-demo kinematics LaTeX table to a file",
    )
    parser.add_argument("--pair-energy-keV", type=float, default=313.0)
    parser.add_argument("--pair-duration-ps", type=float, default=1050.0)
    parser.add_argument("--pair-dt-ps", type=float, default=0.1)
    parser.add_argument(
        "--pair-fine-dt-ps",
        type=float,
        default=None,
        help="optional fine timestep for the initial source-overlap interval",
    )
    parser.add_argument(
        "--pair-fine-duration-ps",
        type=float,
        default=5.0,
        help="duration integrated with --pair-fine-dt-ps before switching to --pair-dt-ps",
    )
    parser.add_argument(
        "--pair-output-dt-ps",
        type=float,
        default=0.1,
        help="trajectory output spacing used by the multirate pair integrator",
    )
    parser.add_argument(
        "--pair-direction",
        type=float,
        nargs=3,
        default=(1.0, 0.0, 0.0),
        metavar=("UX", "UY", "UZ"),
        help="electron launch direction; positron is launched back-to-back",
    )
    parser.add_argument(
        "--pair-origin-m",
        type=float,
        nargs=3,
        default=(0.0, 0.0, 0.0),
        metavar=("X", "Y", "Z"),
    )
    parser.add_argument("--field-span-mm", type=float, default=4.0)
    parser.add_argument("--field-grid", type=int, default=121)
    parser.add_argument("--beambeam-radius-m", type=float, default=0.15)
    parser.add_argument("--beambeam-length-m", type=float, default=0.32)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=DATA_DIR / "boosted_gaussian_witness",
        help="output prefix for pair-demo CSV/PNG files",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.pair_demo:
        return run_pair_demo(args)

    rows = [compute_case(case, args.time_s, args.dt_s) for case in default_cases(args)]

    beta = beta_from_kinetic_energy(args.source_kinetic_MeV)
    gamma = 1.0 / math.sqrt(1.0 - beta * beta)
    print(f"source kinetic energy = {args.source_kinetic_MeV:.6g} MeV")
    print(f"source beta = {beta:.12g}, gamma = {gamma:.12g}")
    print(f"sigma' = {args.sigma_m:.6e} m, Q = {args.charge_C:.6e} C")
    print(f"dt = {args.dt_s:.6e} s")
    print()
    print_text_table(rows)

    table = latex_table(rows)
    if args.latex_table:
        print()
        print(table)
    args.write_latex_table.parent.mkdir(parents=True, exist_ok=True)
    args.write_latex_table.write_text(table + "\n", encoding="utf-8")
    print(f"wrote {args.write_latex_table}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
