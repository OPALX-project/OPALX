#!/usr/bin/env python3
"""Parse and plot the 12-particle BeamBeam trajectory reference file.

The reference file format is documented in TestParticleOrbitSimulation.pptx:

    ----- Test particle Te-1 KIND=2 -----
    t x y s E Px Py Ps Sx Sy Ss

where t, x, y, and s are in meters, E is in eV, P is in eV/c, and S is the
spin vector.  This utility reads the spin columns to preserve the complete
input format, but the current manufactured-solution comparison intentionally
focuses on trajectory and momentum only.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("/Users/adelmann/Desktop/TestParticleOrbit.dat")
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
BOOSTED_GAUSSIAN_PATH = Path(__file__).resolve().parents[1] / "python" / "boosted_gaussian_witness.py"
ELECTRON_REST_EV = 510_998.95
ELEMENTARY_CHARGE_C = 1.602_176_634e-19
DEFAULT_SOURCE_ELECTRONS = 1.25e10
DEFAULT_SOURCE_CHARGE_C = -DEFAULT_SOURCE_ELECTRONS * ELEMENTARY_CHARGE_C
DEFAULT_SOURCE_KINETIC_MEV = 245.0
DEFAULT_MC_SOURCE_PARTICLES = 400_000
DEFAULT_MC_SEED = 20260629
SIGMA_Z_M = 0.6e-3
SIGMA_XY_M = 1.944325075701e-6
C_LIGHT = 299_792_458.0
EPSILON_0 = 8.854_187_812_8e-12

NUMERIC_COLUMNS = ["t", "x", "y", "s", "E", "Px", "Py", "Ps", "Sx", "Sy", "Ss"]
KINEMATIC_COLUMNS = ["t", "x", "y", "s", "E", "Px", "Py", "Ps"]
HEADER_RE = re.compile(
    r"^-+\s+Test particle\s+(?P<name>T(?P<charge_label>e[+-])(?P<pair>\d+))\s+"
    r"KIND=(?P<kind>\d+)\s+-+\s*$"
)
LOG_QUAD_NODES = np.linspace(-36.0, 36.0, 721)


class RigidAnisotropicGaussianSource:
    """Uniformly moving triaxial Gaussian source evaluated in its rest frame."""

    def __init__(
        self,
        name: str,
        charge_C: float,
        sigma_lab_m: tuple[float, float, float],
        beta_z: float,
        center_t0_m: np.ndarray,
        t0_s: float = 0.0,
    ) -> None:
        if abs(beta_z) >= 1.0:
            raise ValueError(f"{name}: |beta_z| must be < 1")
        if any(sigma <= 0.0 for sigma in sigma_lab_m):
            raise ValueError(f"{name}: all source sigmas must be positive")
        self.name = name
        self.charge_C = charge_C
        self.beta = np.array([0.0, 0.0, beta_z], dtype=float)
        self.beta_z = beta_z
        self.gamma = 1.0 / np.sqrt(1.0 - beta_z * beta_z)
        self.center_t0_m = np.asarray(center_t0_m, dtype=float)
        self.t0_s = t0_s

        sigma_x, sigma_y, sigma_z_lab = sigma_lab_m
        # The bunch is specified by lab-frame sigma_z.  In the source rest frame
        # this is dilated by gamma, while transverse sigmas are unchanged.
        self.sigma_rest = np.array([sigma_x, sigma_y, self.gamma * sigma_z_lab], dtype=float)

    def center(self, time_s: float) -> np.ndarray:
        return self.center_t0_m + self.beta * C_LIGHT * (time_s - self.t0_s)


def beta_from_kinetic_energy(kinetic_MeV: float) -> float:
    gamma = 1.0 + kinetic_MeV * 1.0e6 / ELECTRON_REST_EV
    return float(np.sqrt(1.0 - 1.0 / (gamma * gamma)))


def gamma_from_momentum(momentum: np.ndarray) -> float:
    return float(np.sqrt(1.0 + np.dot(momentum, momentum)))


def velocity_from_momentum(momentum: np.ndarray) -> np.ndarray:
    return C_LIGHT * momentum / gamma_from_momentum(momentum)


def source_centroid_offsets(
    sigma_lab_m: tuple[float, float, float],
    n_macroparticles: int,
    seed: int,
) -> dict[str, np.ndarray]:
    """Return reproducible finite-sample source centroid offsets.

    For a Gaussian source represented by N random macroparticles, the sampled
    bunch centroid fluctuates with rms sigma_i / sqrt(N) in each coordinate.
    Offsets for both possible source names are always drawn in a fixed order so
    selecting only one source gives the same offset as selecting it from the
    two-source model.
    """
    if n_macroparticles <= 0:
        zero = np.zeros(3)
        return {"copropagating": zero.copy(), "oncoming": zero.copy()}
    rng = np.random.default_rng(seed)
    rms = np.asarray(sigma_lab_m, dtype=float) / np.sqrt(float(n_macroparticles))
    return {
        "copropagating": rng.normal(0.0, rms),
        "oncoming": rng.normal(0.0, rms),
    }


def anisotropic_gaussian_rest_field(
    r_prime_m: np.ndarray,
    charge_C: float,
    sigma_rest_m: np.ndarray,
) -> np.ndarray:
    """Return E' for a triaxial Gaussian by 1D quadrature.

    The density is normalized as exp(-x_i^2/(2 sigma_i^2)).  The integral is
    the standard ellipsoidal-Gaussian Coulomb field representation.
    """
    if np.allclose(r_prime_m, 0.0, atol=0.0):
        return np.zeros(3)

    sigma2 = sigma_rest_m * sigma_rest_m
    scale = float(np.exp(np.mean(np.log(sigma2))))
    t = scale * np.exp(LOG_QUAD_NODES)
    jacobian = t

    denom_terms = sigma2[None, :] + t[:, None]
    exponent = -0.5 * np.sum((r_prime_m[None, :] * r_prime_m[None, :]) / denom_terms, axis=1)
    common = np.exp(exponent) / np.sqrt(np.prod(denom_terms, axis=1))

    integrals = []
    for axis in range(3):
        integrand = common / denom_terms[:, axis] * jacobian
        integrals.append(float(np.trapezoid(integrand, LOG_QUAD_NODES)))

    coefficient = charge_C / (4.0 * np.pi * EPSILON_0 * np.sqrt(2.0 * np.pi))
    return coefficient * r_prime_m * np.asarray(integrals)


def anisotropic_source_lab_fields(
    position_m: np.ndarray,
    time_s: float,
    source: RigidAnisotropicGaussianSource,
) -> tuple[np.ndarray, np.ndarray]:
    """Return lab-frame E and B from one rigid boosted anisotropic source."""
    displacement = position_m - source.center(time_s)
    r_prime = displacement.copy()
    r_prime[2] *= source.gamma

    e_prime = anisotropic_gaussian_rest_field(r_prime, source.charge_C, source.sigma_rest)
    e_lab = e_prime.copy()
    e_lab[0] *= source.gamma
    e_lab[1] *= source.gamma
    b_lab = np.cross(source.beta * C_LIGHT, e_lab) / (C_LIGHT * C_LIGHT)
    return e_lab, b_lab


def anisotropic_total_lab_fields(
    position_m: np.ndarray,
    time_s: float,
    sources: tuple[RigidAnisotropicGaussianSource, ...],
) -> tuple[np.ndarray, np.ndarray]:
    e_total = np.zeros(3)
    b_total = np.zeros(3)
    for source in sources:
        e_field, b_field = anisotropic_source_lab_fields(position_m, time_s, source)
        e_total += e_field
        b_total += b_field
    return e_total, b_total


def boris_kick(momentum: np.ndarray, e_field: np.ndarray, b_field: np.ndarray, dt_s: float, charge_units: float) -> np.ndarray:
    """Apply OPALX-normalized Boris kick to P=beta*gamma."""
    p = momentum.astype(float).copy()
    p += 0.5 * dt_s * charge_units * C_LIGHT / ELECTRON_REST_EV * e_field

    gamma = gamma_from_momentum(p)
    t_vec = 0.5 * dt_s * charge_units * C_LIGHT * C_LIGHT / (gamma * ELECTRON_REST_EV) * b_field
    w_vec = p + np.cross(p, t_vec)
    s_vec = 2.0 / (1.0 + float(np.dot(t_vec, t_vec))) * t_vec
    p += np.cross(w_vec, s_vec)

    p += 0.5 * dt_s * charge_units * C_LIGHT / ELECTRON_REST_EV * e_field
    return p


def advance_with_fields(
    position: np.ndarray,
    momentum: np.ndarray,
    charge_units: float,
    time_s: float,
    dt_s: float,
    field_function,
    sources,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    e_field, b_field = field_function(position, time_s, sources)
    p_new = boris_kick(momentum, e_field, b_field, dt_s, charge_units)
    v_old = velocity_from_momentum(momentum)
    v_new = velocity_from_momentum(p_new)
    r_new = position + 0.5 * (v_old + v_new) * dt_s
    return r_new, p_new, e_field, b_field


def configure_matplotlib(output_dir: Path) -> None:
    """Configure a writable noninteractive Matplotlib setup before import."""
    cache_root = output_dir / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    cache_root.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)


def parse_reference_file(path: Path) -> pd.DataFrame:
    """Return all trajectory rows from a TestParticleOrbit.dat-style file."""
    rows: list[dict[str, float | int | str]] = []
    current: dict[str, int | str] | None = None
    source_path = path.resolve()
    source_mtime = datetime.fromtimestamp(source_path.stat().st_mtime).astimezone()
    source_mtime_iso = source_mtime.isoformat(timespec="seconds")
    source_date = source_mtime.date().isoformat()

    with source_path.open("r", encoding="utf-8", newline="") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("-----"):
                match = HEADER_RE.match(stripped)
                if match is None:
                    raise ValueError(f"{source_path}:{line_number}: unsupported particle header: {stripped}")
                kind = int(match.group("kind"))
                if kind not in (2, 3):
                    raise ValueError(f"{source_path}:{line_number}: expected KIND=2 or KIND=3, got {kind}")
                current = {
                    "name": match.group("name"),
                    "pair": int(match.group("pair")),
                    "kind": kind,
                    "species": "electron" if kind == 2 else "positron",
                    "charge_units": -1 if kind == 2 else +1,
                    "source_file": str(source_path),
                    "source_mtime": source_mtime_iso,
                    "source_date": source_date,
                }
                continue

            if current is None:
                raise ValueError(f"{source_path}:{line_number}: data row before first particle header")

            values = stripped.split()
            if len(values) != len(NUMERIC_COLUMNS):
                raise ValueError(
                    f"{source_path}:{line_number}: expected {len(NUMERIC_COLUMNS)} columns, got {len(values)}"
                )

            row = dict(current)
            row.update(zip(NUMERIC_COLUMNS, (float(value) for value in values), strict=True))
            rows.append(row)

    if not rows:
        raise ValueError(f"{source_path}: no trajectory rows found")

    data = pd.DataFrame(rows)
    data["step"] = data.groupby("name").cumcount()
    data["t_over_sigma_z"] = data["t"] / SIGMA_Z_M
    data["s_mm"] = 1.0e3 * data["s"]
    data["x_um"] = 1.0e6 * data["x"]
    data["y_um"] = 1.0e6 * data["y"]
    data["kinetic_eV"] = data["E"] - ELECTRON_REST_EV
    data["Px_beta_gamma"] = data["Px"] / ELECTRON_REST_EV
    data["Py_beta_gamma"] = data["Py"] / ELECTRON_REST_EV
    data["Ps_beta_gamma"] = data["Ps"] / ELECTRON_REST_EV
    return data


def particle_summary(data: pd.DataFrame) -> pd.DataFrame:
    """Return one row per particle, excluding spin from the reported metrics."""
    grouped = data.sort_values(["pair", "kind", "step"]).groupby("name", sort=False)
    summary = grouped.agg(
        pair=("pair", "first"),
        kind=("kind", "first"),
        species=("species", "first"),
        charge_units=("charge_units", "first"),
        source_file=("source_file", "first"),
        source_mtime=("source_mtime", "first"),
        source_date=("source_date", "first"),
        rows=("step", "size"),
        t0_m=("t", "first"),
        t1_m=("t", "last"),
        t0_over_sigma_z=("t_over_sigma_z", "first"),
        x0_m=("x", "first"),
        x1_m=("x", "last"),
        y0_m=("y", "first"),
        y1_m=("y", "last"),
        s0_m=("s", "first"),
        s1_m=("s", "last"),
        E0_eV=("E", "first"),
        E1_eV=("E", "last"),
        Px0_eVc=("Px", "first"),
        Px1_eVc=("Px", "last"),
        Py0_eVc=("Py", "first"),
        Py1_eVc=("Py", "last"),
        Ps0_eVc=("Ps", "first"),
        Ps1_eVc=("Ps", "last"),
    )
    return summary.reset_index()


def initial_conditions(data: pd.DataFrame) -> pd.DataFrame:
    """Return the initial state of the 12 reference particles for reuse."""
    first_rows = data.sort_values(["pair", "kind", "step"]).groupby("name", sort=False).head(1)
    columns = [
        "name",
        "pair",
        "kind",
        "species",
        "charge_units",
        "source_file",
        "source_mtime",
        "source_date",
        "t",
        "t_over_sigma_z",
        "x",
        "y",
        "s",
        "E",
        "kinetic_eV",
        "Px",
        "Py",
        "Ps",
        "Px_beta_gamma",
        "Py_beta_gamma",
        "Ps_beta_gamma",
    ]
    return first_rows.loc[:, columns].reset_index(drop=True)


def load_manufactured_module():
    """Load the independent boosted-Gaussian manufactured witness evaluator."""
    spec = importlib.util.spec_from_file_location("boosted_gaussian_witness", BOOSTED_GAUSSIAN_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load manufactured evaluator from {BOOSTED_GAUSSIAN_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_manufactured_sources(
    manufactured,
    source_kinetic_MeV: float,
    source_charge_C: float,
    sigma_m: float,
    source_selection: str,
    source_time_offset_s: float,
) -> tuple[object, object]:
    """Return two equal counter-propagating manufactured electron sources."""
    beta = manufactured.beta_from_kinetic_energy(source_kinetic_MeV)
    source_kwargs = {
        "charge_C": source_charge_C,
        "sigma_m": sigma_m,
        "t0_s": source_time_offset_s,
    }
    copropagating = manufactured.Source(
        "copropagating",
        beta=np.array([0.0, 0.0, +beta]),
        center_t0_m=np.array([0.0, 0.0, 0.0]),
        **source_kwargs,
    )
    oncoming = manufactured.Source(
        "oncoming",
        beta=np.array([0.0, 0.0, -beta]),
        center_t0_m=np.array([0.0, 0.0, 0.0]),
        **source_kwargs,
    )
    if source_selection == "both":
        return (copropagating, oncoming)
    if source_selection == "copropagating":
        return (copropagating,)
    if source_selection == "oncoming":
        return (oncoming,)
    raise ValueError(f"unsupported source selection {source_selection!r}")


def make_anisotropic_sources(
    source_kinetic_MeV: float,
    source_charge_C: float,
    sigma_x_m: float,
    sigma_y_m: float,
    sigma_z_m: float,
    source_selection: str,
    source_time_offset_s: float,
    mc_source_particles: int,
    mc_seed: int,
) -> tuple[RigidAnisotropicGaussianSource, ...]:
    beta = beta_from_kinetic_energy(source_kinetic_MeV)
    centroid_offsets = source_centroid_offsets(
        (sigma_x_m, sigma_y_m, sigma_z_m),
        mc_source_particles,
        mc_seed,
    )
    copropagating = RigidAnisotropicGaussianSource(
        "copropagating",
        charge_C=source_charge_C,
        sigma_lab_m=(sigma_x_m, sigma_y_m, sigma_z_m),
        beta_z=+beta,
        center_t0_m=centroid_offsets["copropagating"],
        t0_s=source_time_offset_s,
    )
    oncoming = RigidAnisotropicGaussianSource(
        "oncoming",
        charge_C=source_charge_C,
        sigma_lab_m=(sigma_x_m, sigma_y_m, sigma_z_m),
        beta_z=-beta,
        center_t0_m=centroid_offsets["oncoming"],
        t0_s=source_time_offset_s,
    )
    if source_selection == "both":
        return (copropagating, oncoming)
    if source_selection == "copropagating":
        return (copropagating,)
    if source_selection == "oncoming":
        return (oncoming,)
    raise ValueError(f"unsupported source selection {source_selection!r}")


def trajectory_error_metrics(reference: pd.DataFrame, manufactured: pd.DataFrame) -> dict[str, float]:
    merged = reference[["name", "step", "x_um"]].merge(
        manufactured[["name", "step", "x_um"]],
        on=["name", "step"],
        suffixes=("_cain", "_manufactured"),
    )
    diff = merged["x_um_manufactured"] - merged["x_um_cain"]
    return {
        "rmse_x_um": float(np.sqrt(np.mean(diff * diff))),
        "mae_x_um": float(np.mean(np.abs(diff))),
        "max_abs_x_um": float(np.max(np.abs(diff))),
    }


def scan_source_phase(reference: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    offsets_ps = np.asarray(args.scan_time_offsets_ps, dtype=float)
    for source_selection in args.scan_source_selection:
        for offset_ps in offsets_ps:
            manufactured = track_manufactured_on_reference_grid(
                reference,
                args.manufactured_model,
                args.source_kinetic_MeV,
                args.source_charge_C,
                args.manufactured_sigma_m,
                args.source_sigma_x_m,
                args.source_sigma_y_m,
                args.source_sigma_z_m,
                args.max_substep_s,
                source_selection,
                offset_ps * 1.0e-12,
                args.mc_source_particles,
                args.mc_seed,
            )
            metrics = trajectory_error_metrics(reference, manufactured)
            rows.append(
                {
                    "manufactured_model": args.manufactured_model,
                    "source_selection": source_selection,
                    "source_time_offset_ps": offset_ps,
                    **metrics,
                }
            )
    return pd.DataFrame(rows).sort_values(["rmse_x_um", "mae_x_um"]).reset_index(drop=True)


def select_best_scan_case(scan: pd.DataFrame) -> tuple[str, float]:
    best = scan.iloc[0]
    return str(best["source_selection"]), float(best["source_time_offset_ps"]) * 1.0e-12


def write_scan_outputs(reference: pd.DataFrame, args: argparse.Namespace) -> None:
    scan = scan_source_phase(reference, args)
    args.scan_csv.parent.mkdir(parents=True, exist_ok=True)
    scan.to_csv(args.scan_csv, index=False)

    best_selection, best_offset_s = select_best_scan_case(scan)
    manufactured = track_manufactured_on_reference_grid(
        reference,
        args.manufactured_model,
        args.source_kinetic_MeV,
        args.source_charge_C,
        args.manufactured_sigma_m,
        args.source_sigma_x_m,
        args.source_sigma_y_m,
        args.source_sigma_z_m,
        args.max_substep_s,
        best_selection,
        best_offset_s,
        args.mc_source_particles,
        args.mc_seed,
    )
    args.scan_comparison_plot.parent.mkdir(parents=True, exist_ok=True)
    plot_cain_and_manufactured(reference, manufactured, args.scan_comparison_plot)
    print("Best scan cases:")
    print(scan.head(10).to_string(index=False))
    print(f"Wrote scan CSV: {args.scan_csv}")
    print(f"Wrote best-scan comparison plot: {args.scan_comparison_plot}")


def track_manufactured_on_reference_grid(
    reference: pd.DataFrame,
    model: str,
    source_kinetic_MeV: float,
    source_charge_C: float,
    sigma_m: float,
    sigma_x_m: float,
    sigma_y_m: float,
    sigma_z_m: float,
    max_substep_s: float,
    source_selection: str,
    source_time_offset_s: float,
    mc_source_particles: int,
    mc_seed: int,
) -> pd.DataFrame:
    """Track the 12 initial particles with the manufactured model on the CAIN ct grid."""
    if max_substep_s <= 0.0:
        raise ValueError("--max-substep-s must be positive")
    if model == "spherical":
        manufactured = load_manufactured_module()
        sources = make_manufactured_sources(
            manufactured,
            source_kinetic_MeV,
            source_charge_C,
            sigma_m,
            source_selection,
            source_time_offset_s,
        )

        def field_function(position_m, time_s, field_sources):
            return manufactured.total_lab_fields(position_m, time_s, field_sources)

        def advance(position, momentum, charge_units, time_s, dt_s):
            return manufactured.advance_witness_boris(position, momentum, charge_units, time_s, dt_s, sources)

        c_light = manufactured.C_LIGHT
        gamma_function = manufactured.gamma_from_momentum
    elif model == "anisotropic":
        sources = make_anisotropic_sources(
            source_kinetic_MeV,
            source_charge_C,
            sigma_x_m,
            sigma_y_m,
            sigma_z_m,
            source_selection,
            source_time_offset_s,
            mc_source_particles,
            mc_seed,
        )

        def field_function(position_m, time_s, field_sources):
            return anisotropic_total_lab_fields(position_m, time_s, field_sources)

        def advance(position, momentum, charge_units, time_s, dt_s):
            return advance_with_fields(
                position,
                momentum,
                charge_units,
                time_s,
                dt_s,
                field_function,
                sources,
            )

        c_light = C_LIGHT
        gamma_function = gamma_from_momentum
    else:
        raise ValueError(f"unsupported manufactured model {model!r}")

    source_offsets = {source.name: source.center_t0_m for source in sources if hasattr(source, "center_t0_m")}
    source_offset_summary = ";".join(
        f"{name}:{offset[0]:.6e},{offset[1]:.6e},{offset[2]:.6e}"
        for name, offset in sorted(source_offsets.items())
    )

    rows: list[dict[str, float | int | str]] = []

    for name, group in reference.sort_values(["pair", "kind", "step"]).groupby("name", sort=False):
        group = group.sort_values("step")
        first = group.iloc[0]
        position = first[["x", "y", "s"]].to_numpy(dtype=float)
        momentum = first[["Px_beta_gamma", "Py_beta_gamma", "Ps_beta_gamma"]].to_numpy(dtype=float)
        charge_units = float(first["charge_units"])

        for local_index, (_row_index, row) in enumerate(group.iterrows()):
            ct_m = float(row["t"])
            time_s = ct_m / c_light
            gamma = gamma_function(momentum)
            e_field, b_field = field_function(position, time_s, sources)
            rows.append(
                {
                    "name": name,
                    "model": model,
                    "source_selection": source_selection,
                    "source_time_offset_s": source_time_offset_s,
                    "mc_source_particles": mc_source_particles,
                    "mc_seed": mc_seed,
                    "source_centroid_offsets_m": source_offset_summary,
                    "pair": int(row["pair"]),
                    "kind": int(row["kind"]),
                    "species": row["species"],
                    "charge_units": int(row["charge_units"]),
                    "source_file": row["source_file"],
                    "source_mtime": row["source_mtime"],
                    "source_date": row["source_date"],
                    "step": int(row["step"]),
                    "t": ct_m,
                    "time_s": time_s,
                    "x": position[0],
                    "y": position[1],
                    "s": position[2],
                    "E": gamma * ELECTRON_REST_EV,
                    "Px": momentum[0] * ELECTRON_REST_EV,
                    "Py": momentum[1] * ELECTRON_REST_EV,
                    "Ps": momentum[2] * ELECTRON_REST_EV,
                    "t_over_sigma_z": ct_m / SIGMA_Z_M,
                    "s_mm": 1.0e3 * position[2],
                    "x_um": 1.0e6 * position[0],
                    "y_um": 1.0e6 * position[1],
                    "kinetic_eV": (gamma - 1.0) * ELECTRON_REST_EV,
                    "Px_beta_gamma": momentum[0],
                    "Py_beta_gamma": momentum[1],
                    "Ps_beta_gamma": momentum[2],
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
            if local_index == len(group) - 1:
                continue
            next_ct_m = float(group.iloc[local_index + 1]["t"])
            interval_s = (next_ct_m - ct_m) / c_light
            if interval_s < 0.0:
                raise ValueError(f"{name}: reference ct grid is not monotonic at step {row['step']}")
            n_substeps = max(1, int(np.ceil(interval_s / max_substep_s)))
            substep_s = interval_s / n_substeps
            for substep in range(n_substeps):
                position, momentum, _e_field, _b_field = advance(
                    position, momentum, charge_units, time_s + substep * substep_s, substep_s
                )

    manufactured_data = pd.DataFrame(rows)
    return manufactured_data


def plot_x_vs_s(data: pd.DataFrame, output: Path) -> None:
    """Reproduce the slide's x-versus-s trajectory plot."""
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 13,
            "axes.labelsize": 22,
            "axes.titlesize": 21,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
        }
    )

    line_styles = [
        "-",
        (0, (8, 6)),
        (0, (5, 5)),
        (0, (10, 8)),
        (0, (2, 3)),
        (0, (9, 4, 2, 4)),
    ]

    fig, ax = plt.subplots(figsize=(8.0, 6.0), dpi=180)
    for (_name, pair, kind), group in data.groupby(["name", "pair", "kind"], sort=False):
        color = "#c00000" if kind == 2 else "#0017b8"
        ax.plot(
            group["s_mm"],
            group["x_um"],
            color=color,
            linewidth=0.9,
            linestyle=line_styles[(int(pair) - 1) % len(line_styles)],
        )

    ax.axhline(0.0, color="black", linewidth=0.8, linestyle=(0, (2, 6)))
    ax.set_title("Test Particle Trajectory")
    ax.set_xlabel("S (mm)")
    ax.set_ylabel(r"x ($\mu$m)")
    ax.set_xlim(-0.75, 1.65)
    ax.set_ylim(-9.0, 20.0)
    ax.tick_params(direction="out", top=True, right=True, which="both")
    ax.minorticks_on()
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_cain_and_manufactured(reference: pd.DataFrame, manufactured: pd.DataFrame, output: Path) -> None:
    """Plot CAIN/reference and manufactured x(s) trajectories side by side."""
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    line_styles = [
        "-",
        (0, (8, 6)),
        (0, (5, 5)),
        (0, (10, 8)),
        (0, (2, 3)),
        (0, (9, 4, 2, 4)),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6), dpi=180, sharex=True, sharey=True)
    model_label = str(manufactured["model"].iloc[0]) if "model" in manufactured and not manufactured.empty else "manufactured"
    panels = (
        (axes[0], reference, "CAIN reference"),
        (axes[1], manufactured, f"Manufactured tracking ({model_label})"),
    )
    for ax, frame, title in panels:
        for (_name, pair, kind), group in frame.groupby(["name", "pair", "kind"], sort=False):
            color = "#c00000" if kind == 2 else "#0017b8"
            ax.plot(
                group["s_mm"],
                group["x_um"],
                color=color,
                linewidth=0.9,
                linestyle=line_styles[(int(pair) - 1) % len(line_styles)],
            )
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle=(0, (2, 6)))
        ax.set_title(title)
        ax.set_xlabel("S (mm)")
        ax.set_xlim(-0.75, 1.65)
        ax.set_ylim(-9.0, 20.0)
        ax.tick_params(direction="out", top=True, right=True, which="both")
        ax.minorticks_on()

    axes[0].set_ylabel(r"x ($\mu$m)")
    fig.suptitle("12 Test Particle Trajectories")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_coordinate_comparison(
    reference: pd.DataFrame,
    manufactured: pd.DataFrame,
    output: Path,
    coordinate: str,
) -> None:
    """Plot CAIN/reference and manufactured coordinate trajectories side by side."""
    if coordinate == "y":
        x_column = "s_mm"
        y_column = "y_um"
        xlabel = "S (mm)"
        ylabel = r"y ($\mu$m)"
        suptitle = "12 Test Particle y Trajectories"
        xlim = (-0.75, 1.65)
        ylim = None
    elif coordinate == "z":
        x_column = "t"
        y_column = "s_mm"
        xlabel = "ct (mm)"
        ylabel = "z = S (mm)"
        suptitle = "12 Test Particle z Trajectories"
        xlim = (-0.95, 1.85)
        ylim = (-0.95, 1.75)
    else:
        raise ValueError(f"unsupported coordinate {coordinate!r}")

    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 12,
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    line_styles = [
        "-",
        (0, (8, 6)),
        (0, (5, 5)),
        (0, (10, 8)),
        (0, (2, 3)),
        (0, (9, 4, 2, 4)),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6), dpi=180, sharex=True, sharey=True)
    model_label = str(manufactured["model"].iloc[0]) if "model" in manufactured and not manufactured.empty else "manufactured"
    panels = (
        (axes[0], reference, "CAIN reference"),
        (axes[1], manufactured, f"Manufactured tracking ({model_label})"),
    )
    for ax, frame, title in panels:
        for (_name, pair, kind), group in frame.groupby(["name", "pair", "kind"], sort=False):
            color = "#c00000" if kind == 2 else "#0017b8"
            x_values = 1.0e3 * group[x_column] if coordinate == "z" else group[x_column]
            ax.plot(
                x_values,
                group[y_column],
                color=color,
                linewidth=0.9,
                linestyle=line_styles[(int(pair) - 1) % len(line_styles)],
            )
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle=(0, (2, 6)))
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.tick_params(direction="out", top=True, right=True, which="both")
        ax.minorticks_on()

    axes[0].set_ylabel(ylabel)
    fig.suptitle(suptitle)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_manufactured_x_vs_s(manufactured: pd.DataFrame, output: Path) -> None:
    """Plot only the manufactured x-versus-s trajectories."""
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 13,
            "axes.labelsize": 22,
            "axes.titlesize": 21,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
        }
    )
    line_styles = [
        "-",
        (0, (8, 6)),
        (0, (5, 5)),
        (0, (10, 8)),
        (0, (2, 3)),
        (0, (9, 4, 2, 4)),
    ]

    fig, ax = plt.subplots(figsize=(8.0, 6.0), dpi=180)
    for (_name, pair, kind), group in manufactured.groupby(["name", "pair", "kind"], sort=False):
        color = "#c00000" if kind == 2 else "#0017b8"
        ax.plot(
            group["s_mm"],
            group["x_um"],
            color=color,
            linewidth=0.9,
            linestyle=line_styles[(int(pair) - 1) % len(line_styles)],
        )

    ax.axhline(0.0, color="black", linewidth=0.8, linestyle=(0, (2, 6)))
    model_label = str(manufactured["model"].iloc[0]) if "model" in manufactured and not manufactured.empty else "manufactured"
    ax.set_title(f"Manufactured Test Particle Trajectories ({model_label})")
    ax.set_xlabel("S (mm)")
    ax.set_ylabel(r"x ($\mu$m)")
    ax.set_xlim(-0.75, 1.65)
    ax.set_ylim(-9.0, 20.0)
    ax.tick_params(direction="out", top=True, right=True, which="both")
    ax.minorticks_on()
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def write_outputs(data: pd.DataFrame, args: argparse.Namespace) -> None:
    """Write requested tables and plots."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manufactured_data = None

    if args.trajectory_csv is not None:
        trajectory_columns = [
            "name",
            "pair",
            "kind",
            "species",
            "charge_units",
            "source_file",
            "source_mtime",
            "source_date",
            "step",
            *KINEMATIC_COLUMNS,
            "t_over_sigma_z",
            "s_mm",
            "x_um",
            "y_um",
            "kinetic_eV",
            "Px_beta_gamma",
            "Py_beta_gamma",
            "Ps_beta_gamma",
        ]
        args.trajectory_csv.parent.mkdir(parents=True, exist_ok=True)
        data.loc[:, trajectory_columns].to_csv(args.trajectory_csv, index=False)

    if args.initial_csv is not None:
        args.initial_csv.parent.mkdir(parents=True, exist_ok=True)
        initial_conditions(data).to_csv(args.initial_csv, index=False)

    if args.summary_csv is not None:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        particle_summary(data).to_csv(args.summary_csv, index=False)

    if args.plot is not None:
        plot_x_vs_s(data, args.plot)

    if (
        args.manufactured_csv is not None
        or args.comparison_plot is not None
        or args.comparison_y_plot is not None
        or args.comparison_z_plot is not None
        or args.manufactured_plot is not None
    ):
        manufactured_data = track_manufactured_on_reference_grid(
            data,
            args.manufactured_model,
            args.source_kinetic_MeV,
            args.source_charge_C,
            args.manufactured_sigma_m,
            args.source_sigma_x_m,
            args.source_sigma_y_m,
            args.source_sigma_z_m,
            args.max_substep_s,
            args.source_selection,
            args.source_time_offset_s,
            args.mc_source_particles,
            args.mc_seed,
        )

    if args.manufactured_csv is not None:
        args.manufactured_csv.parent.mkdir(parents=True, exist_ok=True)
        manufactured_data.to_csv(args.manufactured_csv, index=False)

    if args.comparison_plot is not None:
        plot_cain_and_manufactured(data, manufactured_data, args.comparison_plot)

    if args.comparison_y_plot is not None:
        plot_coordinate_comparison(data, manufactured_data, args.comparison_y_plot, "y")

    if args.comparison_z_plot is not None:
        plot_coordinate_comparison(data, manufactured_data, args.comparison_z_plot, "z")

    if args.manufactured_plot is not None:
        plot_manufactured_x_vs_s(manufactured_data, args.manufactured_plot)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"TestParticleOrbit.dat path. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--trajectory-csv",
        type=Path,
        default=None,
        help="Full trajectory CSV path. Defaults to OUTPUT_DIR/track12_trajectory.csv.",
    )
    parser.add_argument(
        "--initial-csv",
        type=Path,
        default=None,
        help="Initial-condition CSV path. Defaults to OUTPUT_DIR/track12_initial_conditions.csv.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=None,
        help="Particle summary CSV path. Defaults to OUTPUT_DIR/track12_summary.csv.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="x(s) trajectory plot path. Defaults to OUTPUT_DIR/track12_x_vs_s.png.",
    )
    parser.add_argument(
        "--manufactured-csv",
        type=Path,
        default=None,
        help="Manufactured trajectory CSV path. Defaults to OUTPUT_DIR/track12_manufactured_trajectory.csv.",
    )
    parser.add_argument(
        "--comparison-plot",
        type=Path,
        default=None,
        help="CAIN/reference versus manufactured x(s) plot. Defaults to OUTPUT_DIR/track12_cain_vs_manufactured.png.",
    )
    parser.add_argument(
        "--comparison-y-plot",
        type=Path,
        default=None,
        help="CAIN/reference versus manufactured y(s) plot.",
    )
    parser.add_argument(
        "--comparison-z-plot",
        type=Path,
        default=None,
        help="CAIN/reference versus manufactured z(ct) plot.",
    )
    parser.add_argument(
        "--manufactured-plot",
        type=Path,
        default=None,
        help="Manufactured-only x(s) plot. Defaults to OUTPUT_DIR/track12_manufactured_x_vs_s.png.",
    )
    parser.add_argument(
        "--source-kinetic-MeV",
        type=float,
        default=DEFAULT_SOURCE_KINETIC_MEV,
        help="Manufactured source electron kinetic energy.",
    )
    parser.add_argument(
        "--source-charge-C",
        type=float,
        default=DEFAULT_SOURCE_CHARGE_C,
        help="Manufactured source bunch charge in C.",
    )
    parser.add_argument(
        "--manufactured-sigma-m",
        type=float,
        default=SIGMA_Z_M,
        help="Manufactured spherical Gaussian rest-frame sigma in m.",
    )
    parser.add_argument(
        "--manufactured-model",
        choices=("anisotropic", "spherical"),
        default="anisotropic",
        help="Manufactured source model. Default: anisotropic.",
    )
    parser.add_argument(
        "--source-sigma-x-m",
        type=float,
        default=SIGMA_XY_M,
        help="Lab-frame anisotropic source sigma_x in m.",
    )
    parser.add_argument(
        "--source-sigma-y-m",
        type=float,
        default=SIGMA_XY_M,
        help="Lab-frame anisotropic source sigma_y in m.",
    )
    parser.add_argument(
        "--source-sigma-z-m",
        type=float,
        default=SIGMA_Z_M,
        help="Lab-frame anisotropic source sigma_z in m.",
    )
    parser.add_argument(
        "--max-substep-s",
        type=float,
        default=1.0e-15,
        help="Maximum manufactured tracker timestep between CAIN output samples.",
    )
    parser.add_argument(
        "--source-selection",
        choices=("both", "oncoming", "copropagating"),
        default="oncoming",
        help="Source bunches included in the manufactured force model.",
    )
    parser.add_argument(
        "--source-time-offset-s",
        type=float,
        default=0.0,
        help="Time at which source centers are at the IP in the manufactured model.",
    )
    parser.add_argument(
        "--mc-source-particles",
        type=int,
        default=DEFAULT_MC_SOURCE_PARTICLES,
        help="Macroparticle count used for 1/sqrt(N) source centroid noise; use 0 to disable.",
    )
    parser.add_argument(
        "--mc-seed",
        type=int,
        default=DEFAULT_MC_SEED,
        help="Random seed for reproducible source centroid noise.",
    )
    parser.add_argument(
        "--scan-source-phase",
        action="store_true",
        help="Scan source selection and source time offsets, then plot the best case.",
    )
    parser.add_argument(
        "--scan-source-selection",
        nargs="+",
        choices=("both", "oncoming", "copropagating"),
        default=("both", "oncoming"),
        help="Source selections used by --scan-source-phase.",
    )
    parser.add_argument(
        "--scan-time-offsets-ps",
        type=float,
        nargs="+",
        default=(-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0),
        help="Source time offsets in ps used by --scan-source-phase.",
    )
    parser.add_argument(
        "--scan-csv",
        type=Path,
        default=None,
        help="Scan metrics CSV path. Defaults to OUTPUT_DIR/track12_source_phase_scan.csv.",
    )
    parser.add_argument(
        "--scan-comparison-plot",
        type=Path,
        default=None,
        help="Best scan CAIN/manufactured comparison plot.",
    )
    parser.add_argument("--no-csv", action="store_true", help="Do not write CSV outputs.")
    parser.add_argument("--no-plot", action="store_true", help="Do not write the x(s) plot.")
    parser.add_argument(
        "--no-manufactured",
        action="store_true",
        help="Do not write manufactured trajectory or CAIN-versus-manufactured plot.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir = args.output_dir.resolve()

    if args.trajectory_csv is None and not args.no_csv:
        args.trajectory_csv = args.output_dir / "track12_trajectory.csv"
    if args.initial_csv is None and not args.no_csv:
        args.initial_csv = args.output_dir / "track12_initial_conditions.csv"
    if args.summary_csv is None and not args.no_csv:
        args.summary_csv = args.output_dir / "track12_summary.csv"
    if args.plot is None and not args.no_plot:
        args.plot = args.output_dir / "track12_x_vs_s.png"
    if args.manufactured_csv is None and not args.no_csv and not args.no_manufactured:
        args.manufactured_csv = args.output_dir / "track12_manufactured_trajectory.csv"
    if args.comparison_plot is None and not args.no_plot and not args.no_manufactured:
        args.comparison_plot = args.output_dir / "track12_cain_vs_manufactured.png"
    if args.comparison_y_plot is None and not args.no_plot and not args.no_manufactured:
        args.comparison_y_plot = args.output_dir / "track12_cain_vs_manufactured_y.png"
    if args.comparison_z_plot is None and not args.no_plot and not args.no_manufactured:
        args.comparison_z_plot = args.output_dir / "track12_cain_vs_manufactured_z.png"
    if args.manufactured_plot is None and not args.no_plot and not args.no_manufactured:
        args.manufactured_plot = args.output_dir / "track12_manufactured_x_vs_s.png"
    if args.scan_csv is None:
        args.scan_csv = args.output_dir / "track12_source_phase_scan.csv"
    if args.scan_comparison_plot is None:
        args.scan_comparison_plot = args.output_dir / "track12_source_phase_scan_best.png"

    if args.no_csv:
        args.trajectory_csv = None
        args.initial_csv = None
        args.summary_csv = None
        args.manufactured_csv = None
    if args.no_plot:
        args.plot = None
        args.comparison_plot = None
        args.comparison_y_plot = None
        args.comparison_z_plot = None
        args.manufactured_plot = None
    if args.no_manufactured:
        args.manufactured_csv = None
        args.comparison_plot = None
        args.comparison_y_plot = None
        args.comparison_z_plot = None
        args.manufactured_plot = None

    data = parse_reference_file(args.input)
    summary = particle_summary(data)
    write_outputs(data, args)
    if args.scan_source_phase:
        write_scan_outputs(data, args)

    source_mtime = data["source_mtime"].iloc[0]
    source_date = data["source_date"].iloc[0]
    print(
        f"Read {len(data)} rows for {summary.shape[0]} particles from {args.input} "
        f"(source_date={source_date}, source_mtime={source_mtime})"
    )
    print(summary.loc[:, ["name", "kind", "species", "rows", "t0_over_sigma_z", "x0_m", "s0_m", "s1_m"]].to_string(index=False))
    if args.trajectory_csv is not None:
        print(f"Wrote trajectory CSV: {args.trajectory_csv}")
    if args.initial_csv is not None:
        print(f"Wrote initial-condition CSV: {args.initial_csv}")
    if args.summary_csv is not None:
        print(f"Wrote summary CSV: {args.summary_csv}")
    if args.plot is not None:
        print(f"Wrote plot: {args.plot}")
    if args.manufactured_csv is not None:
        print(f"Wrote manufactured trajectory CSV: {args.manufactured_csv}")
    if args.comparison_plot is not None:
        print(f"Wrote CAIN/manufactured comparison plot: {args.comparison_plot}")
    if args.comparison_y_plot is not None:
        print(f"Wrote CAIN/manufactured y comparison plot: {args.comparison_y_plot}")
    if args.comparison_z_plot is not None:
        print(f"Wrote CAIN/manufactured z comparison plot: {args.comparison_z_plot}")
    if args.manufactured_plot is not None:
        print(f"Wrote manufactured plot: {args.manufactured_plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
