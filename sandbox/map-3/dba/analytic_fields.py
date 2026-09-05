#!/usr/bin/env python3
"""Reproduce native SBEND fields on the nominal map-3 design line.

This is an analytic field model, not a tracked orbit or an exact transfer map.
Source authority: src/AbsBeamline/BendFieldModel.h and SBend::makeFieldInputs.
Detailed equations and limitations: physics/map-3-fringe-dba.qmd in the manual.
All field vectors use the local design-line (x, y, s) basis, in tesla.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
C_LIGHT = 299792458.0
ENGE = np.array([0.478959, 1.911289, -1.185953, 1.630554, -1.082657, 0.318111])


@dataclass(frozen=True)
class Parameters:
    momentum_GeV_c: float
    charge_e: float
    radius_m: float
    body_angle_deg: float
    hgap_m: float
    fint: float
    field_free_drift_m: float
    quadrupole_length_m: float
    quadrupole_k1_seed_m2: float
    probe_y_m: float

    def __post_init__(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError("All parameters must be finite")
        if min(self.momentum_GeV_c, self.radius_m, self.body_angle_deg,
               self.field_free_drift_m, self.quadrupole_length_m) <= 0:
            raise ValueError("Momentum, lengths and angle must be positive")
        if self.charge_e == 0 or self.hgap_m < 0:
            raise ValueError("Charge must be nonzero and HGAP nonnegative")
        if self.hgap_m > 0 and abs(self.probe_y_m) >= self.hgap_m:
            raise ValueError("The diagnostic probe must be inside the vertical aperture")

    @property
    def body_length(self):
        return self.radius_m * np.deg2rad(self.body_angle_deg)

    @property
    def gap(self):
        return 2 * self.hgap_m

    @property
    def half_width(self):
        return 5 * self.gap

    @property
    def b0(self):
        return self.momentum_GeV_c * 1e9 / (self.charge_e * C_LIGHT * self.radius_m)

    @property
    def quad_gradient(self):
        # Standalone OpalQuadrupole/Multipole does NOT divide by charge,
        # unlike SBend::makeFieldInputs. Preserve the current source convention.
        return self.momentum_GeV_c * 1e9 / C_LIGHT * self.quadrupole_k1_seed_m2


def load_parameters(path=ROOT / "parameters.json"):
    return Parameters(**json.loads(Path(path).read_text()))


def layout(p):
    """Nonoverlapping half-open support intervals on the nominal design line."""
    rows = []
    start = 0.0
    for name, kind, length in (
        ("B1", "bend", p.body_length + 2 * p.half_width),
        ("D1", "drift", p.field_free_drift_m),
        ("QACH", "quadrupole", p.quadrupole_length_m),
        ("D2", "drift", p.field_free_drift_m),
        ("B2", "bend", p.body_length + 2 * p.half_width),
    ):
        body_start = start + (p.half_width if kind == "bend" else 0)
        body_length = p.body_length if kind == "bend" else length
        rows.append(dict(name=name, kind=kind, support_start_m=start,
                         support_end_m=start + length, body_start_m=body_start,
                         body_end_m=body_start + body_length))
        start += length
    return pd.DataFrame(rows)


def enge(d, gap):
    """Return F(d), dF/dd, d²F/dd²; same saturation thresholds as the C++."""
    d = np.asarray(d, dtype=float)
    if gap <= 0:
        return np.ones_like(d), np.zeros_like(d), np.zeros_like(d)
    u = d / gap
    polynomial = np.polynomial.polynomial.polyval(u, ENGE)
    first = np.polynomial.polynomial.polyval(u, np.arange(1, 6) * ENGE[1:]) / gap
    second = np.polynomial.polynomial.polyval(
        u, np.arange(2, 6) * np.arange(1, 5) * ENGE[2:]) / gap**2
    exponential = np.exp(np.clip(polynomial, -80, 80))
    f = 1 / (1 + exponential)
    # Use the exponential form to retain small derivatives where F rounds to 1.
    fp = -first * exponential * f**2
    fpp = -(second + first**2) * exponential * f**2 + 2 * first**2 * exponential**2 * f**3
    saturated = abs(polynomial) > 80
    f = np.where(polynomial > 80, 0, np.where(polynomial < -80, 1, f))
    return f, np.where(saturated, 0, fp), np.where(saturated, 0, fpp)


def envelope(z, length, gap):
    """Return A, dA/dz, d²A/dz², and active face (entrance=0, exit=1)."""
    entrance, exit_ = enge(-np.asarray(z), gap), enge(np.asarray(z) - length, gap)
    use_exit = exit_[0] < entrance[0]  # A tie selects entrance, as in C++.
    return (np.where(use_exit, exit_[0], entrance[0]),
            np.where(use_exit, exit_[1], -entrance[1]),
            np.where(use_exit, exit_[2], entrance[2]), use_exit.astype(int))


def bend_field(z, x, y, p, fint=None):
    """Pure sector bend, including external support mask and distributed FINT."""
    z, x, y = np.broadcast_arrays(np.asarray(z, float), x, y)
    active = (z >= -p.half_width) & (z < p.body_length + p.half_width)
    # Evaluate only inside the support; no unnecessary polynomial extrapolation.
    zs = np.clip(z, -p.half_width, p.body_length + p.half_width)
    a, ap, app, face = envelope(zs, p.body_length, p.gap)
    coefficient = 0.0
    if p.gap > 0:
        h = 1 / p.radius_m
        psi = h * p.hgap_m * (p.fint if fint is None else fint)
        ky = h * np.tan(psi)  # E1=E2=0; no pole-face rotations.
        span = abs(enge(-p.half_width, p.gap)[0] - enge(p.half_width, p.gap)[0])
        coefficient = p.b0 / h * ky / span
    field = np.stack((np.where(face == 0, coefficient, -coefficient) * ap * y,
                      p.b0 * (a - 0.5 * app * y**2), p.b0 * ap * y), axis=-1)
    return np.where(active[..., None], field, 0.0)


def lattice_field(s, x, y, p, fint=None):
    """Sample a fixed local offset, not a particle trajectory; drift fields zero."""
    s, x, y = np.broadcast_arrays(np.asarray(s, float), x, y)
    result = np.zeros(s.shape + (3,))
    for row in layout(p).itertuples():
        if row.kind == "bend":
            result += bend_field(s - row.body_start_m, x, y, p, fint)
        elif row.kind == "quadrupole":
            active = (s >= row.body_start_m) & (s < row.body_end_m)
            result[..., 0] += active * p.quad_gradient * y
            result[..., 1] += active * p.quad_gradient * x
    return result


def sample_fields(p):
    spans = layout(p)
    boundaries = spans[["support_start_m", "support_end_m", "body_start_m", "body_end_m"]].to_numpy().ravel()
    length = spans.support_end_m.iloc[-1]
    s = np.unique(np.concatenate((np.linspace(0, length, 6001), boundaries,
                                  np.nextafter(boundaries, -np.inf),
                                  np.nextafter(boundaries, np.inf))))
    s = s[(s >= 0) & (s <= length)]
    centre = lattice_field(s, 0, 0, p)
    probe = lattice_field(s, 0, p.probe_y_m, p)
    no_fint = lattice_field(s, 0, p.probe_y_m, p, fint=0)
    frame = pd.DataFrame({"s_m": s})
    for i, component in enumerate(("Bx", "By", "Bs")):
        frame[f"{component}_axis_T"] = centre[:, i]
        frame[f"{component}_probe_T"] = probe[:, i]
        frame[f"{component}_probe_fint0_T"] = no_fint[:, i]
    frame["delta_Bx_fint_T"] = probe[:, 0] - no_fint[:, 0]
    return frame


def draw_layout(ax, spans):
    for row in spans.itertuples():
        color = {"bend": "#245b8a", "drift": "#6e7378", "quadrupole": "#a24727"}[row.kind]
        a, b = row.support_start_m, row.support_end_m
        ax.hlines(0.55, a, b, color=color, lw=1.4)
        ax.vlines([a, b], 0.43, 0.67, color=color, lw=1.4)
        ax.text((a + b) / 2, 0.85, row.name, ha="center", va="bottom", color=color)
        if row.kind == "bend":
            ax.hlines(0.22, row.body_start_m, row.body_end_m, color=color, lw=5)
            ax.text((a + b) / 2, -0.02, "body", ha="center", va="top", fontsize=10)
    ax.set_ylim(-0.4, 1.4)
    ax.axis("off")


def decorate(ax, spans):
    for row in spans.itertuples():
        if row.kind == "bend":
            ax.axvspan(row.body_start_m, row.body_end_m, color="#245b8a", alpha=0.07, lw=0)
        ax.axvline(row.support_start_m, color="0.8", lw=0.65, zorder=0)
    ax.grid(axis="y", color="0.9", lw=0.6)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_xlim(0, spans.support_end_m.iloc[-1])


def plot_fields(data, p, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                         "axes.labelsize": 13, "axes.titlesize": 13,
                         "savefig.dpi": 220, "svg.fonttype": "none"})
    spans = layout(p)
    s = data.s_m.to_numpy()
    fig = plt.figure(figsize=(11.5, 5.4), layout="constrained")
    axes = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [1, 4]})
    draw_layout(axes[0], spans)
    ax = axes[1]
    ax.plot(s, data.By_axis_T, color="#245b8a", lw=1.9, label=r"$B_y$")
    ax.plot(s, data.Bx_axis_T, color="#ba5b2a", lw=1.6, ls="--", label=r"$B_x=0$")
    ax.plot(s, data.Bs_axis_T, color="#27806c", lw=1.4, ls=":", label=r"$B_s=0$")
    decorate(ax, spans)
    ax.set_ylabel("Magnetic field [T]")
    ax.set_xlabel(r"Nominal design-line coordinate $s$ [m]")
    ax.legend(loc="lower center", ncol=3, frameon=False)
    ax.set_ylim(min(0, p.b0) - 0.15 * abs(p.b0), max(0, p.b0) + 0.15 * abs(p.b0))
    fig.suptitle("map-3 DBA — analytic centreline field\n"
                 f"HGAP = {p.hgap_m:g} m, FINT = {p.fint:g}; support brackets and shaded bend bodies",
                 fontsize=14)
    for suffix in ("png", "svg"):
        fig.savefig(output / f"magnetic-field-vs-s.{suffix}")
    plt.close(fig)

    fig = plt.figure(figsize=(11.5, 9.0), layout="constrained")
    axes = fig.subplots(4, 1, sharex=True, gridspec_kw={"height_ratios": [0.65, 2, 2, 2]})
    draw_layout(axes[0], spans)
    for ax in axes[1:]:
        decorate(ax, spans)
    axes[1].plot(s, data.Bx_probe_T * 1e3, color="#ba5b2a", lw=1.7, label=f"FINT = {p.fint:g}")
    axes[1].plot(s, data.Bx_probe_fint0_T * 1e3, color="0.3", ls="--", lw=1, label="FINT = 0")
    axes[1].set_ylabel(r"$B_x$ [mT]")
    axes[1].set_title("Full radial field (quadrupole strength is an unmatched seed)", loc="left")
    axes[1].legend(loc="lower left", ncol=2, frameon=False)
    axes[2].plot(s, data.delta_Bx_fint_T * 1e6, color="#ba5b2a", lw=1.7)
    axes[2].set_ylabel(r"$\Delta B_x$ [$\mu$T]")
    axes[2].set_title(rf"FINT contribution only: $B_x({p.fint:g})-B_x(0)$", loc="left")
    axes[3].plot(s, data.Bs_probe_T * 1e3, color="#27806c", lw=1.7)
    axes[3].set_ylabel(r"$B_s$ [mT]")
    axes[3].set_title("Longitudinal fringe field (unchanged by FINT)", loc="left")
    axes[3].set_xlabel(r"Nominal design-line coordinate $s$ [m]")
    fig.suptitle(f"map-3 DBA — fixed off-axis probe, x = 0, y = {p.probe_y_m * 1e3:g} mm\n"
                 "Local (x, y, s) field components; not a tracked orbit", fontsize=14)
    for suffix in ("png", "svg"):
        fig.savefig(output / f"off-axis-fringe-fields.{suffix}")
    plt.close(fig)


def summary(p, parameter_path):
    effective, error = quad(lambda z: float(envelope(z, p.body_length, p.gap)[0]),
                            -p.half_width, p.body_length + p.half_width,
                            points=[0, p.body_length / 2, p.body_length], epsabs=1e-12, epsrel=1e-12)
    sources = [Path(__file__), parameter_path, REPO / "src/AbsBeamline/BendFieldModel.h",
               REPO / "src/AbsBeamline/SBend.cpp", REPO / "src/Elements/OpalSBend.cpp",
               REPO / "src/Elements/OpalQuadrupole.cpp", REPO / "src/AbsBeamline/Multipole.cpp",
               REPO / "src/Elements/PlacementResolver.cpp", REPO / "src/BeamlineGeometry/Geometry.h"]
    return dict(model="native SBEND field law on nominal design line; no tracked orbit",
                matched_achromat=False, simulation_run=False, parameters=asdict(p),
                full_gap_m=p.gap, fringe_half_width_m=p.half_width,
                total_length_m=float(layout(p).support_end_m.iloc[-1]),
                body_length_m=p.body_length, dipole_field_T=p.b0,
                quadrupole_gradient_seed_T_m=p.quad_gradient,
                envelope_integral_m=effective, quadrature_error_estimate_m=error,
                integrated_By_T_m=p.b0 * effective,
                nominal_line_curvature_integral_deg=float(np.rad2deg(effective / p.radius_m)),
                note="Curvature integral is NOT the actual orbit deflection. No dipole renormalization.",
                sources_sha256={str(f.relative_to(REPO)) if f.is_relative_to(REPO) else str(f):
                                hashlib.sha256(f.read_bytes()).hexdigest() for f in sources})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=ROOT / "parameters.json")
    parser.add_argument("--output", type=Path, default=ROOT / "analytic")
    args = parser.parse_args()
    p = load_parameters(args.parameters)
    args.output.mkdir(parents=True, exist_ok=True)
    data = sample_fields(p)
    data.to_csv(args.output / "fields-vs-s.csv", index=False, float_format="%.15e")
    spans = layout(p)
    spans.to_csv(args.output / "layout.csv", index=False, float_format="%.15e")
    report = summary(p, args.parameters.resolve())
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_fields(data, p, args.output)
    print(spans.to_string(index=False, float_format=lambda x: f"{x:.9f}"))
    print(json.dumps({key: value for key, value in report.items() if key != "sources_sha256"}, indent=2))


if __name__ == "__main__":
    main()
