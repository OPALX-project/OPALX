#!/usr/bin/env python3
"""Independent linear map of the analytic map-3 field model.

Integrate analytic Lorentz sensitivities, not finite-difference shadow rays.
Path length is the independent variable; quadrupole face crossing sensitivity
uses saltation, and the final map is projected onto the reference-normal plane.
This is a numerically evaluated analytic-field reference, NOT a closed form.
Only the horizontal mid-plane reference is supported by this field Jacobian.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-map3-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/opalx-map3-cache")
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import solve_ivp

from analytic_fields import C_LIGHT, ENGE, ROOT, REPO, layout, load_parameters
from validate_opalx_fields import placed_elements, rotation_y


def cross_matrix(v):
    x, y, z = v
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def particle_constants():
    """Read the actual input momentum and source rest mass, independently of P0."""
    parts = (ROOT.parent / "reference-particle.txt").read_text().splitlines()[2].split()
    assert [float(parts[i]) for i in range(5)] == [0.0] * 5
    u0 = float(parts[5])
    text = (REPO / "src/Physics/Physics.h").read_text()
    mass = float(re.search(r"constexpr double m_e\s*=\s*([^;]+);", text).group(1))
    return u0, mass


def enge_first(d, gap):
    """Scalar Enge function and analytic first derivative in metres."""
    if gap == 0:
        return 1.0, 0.0
    u = d / gap
    exponent = float(ENGE[5])
    for coefficient in ENGE[4::-1]:
        exponent = exponent * u + float(coefficient)
    if exponent > 80:
        return 0.0, 0.0
    if exponent < -80:
        return 1.0, 0.0
    prime = 5 * float(ENGE[5])
    for power in range(4, 0, -1):
        prime = prime * u + power * float(ENGE[power])
    e = math.exp(exponent)
    f = 1 / (1 + e)
    return f, -(prime / gap) * e * f * f


@dataclass
class Face:
    name: str
    origin: np.ndarray
    normal: np.ndarray
    next_quad_on: bool


class PlanarField:
    """Analytic B and dB/dr on y=0, with a prescribed quad mode per interval."""
    def __init__(self, parameters):
        self.p = parameters
        self.poses = list(placed_elements(parameters))
        self.theta = math.radians(parameters.body_angle_deg)
        self.sine, self.cosine = math.sin(self.theta), math.cos(self.theta)
        self.edge_c = 0.0
        if parameters.gap:
            span = abs(enge_first(-parameters.half_width, parameters.gap)[0]
                       - enge_first(parameters.half_width, parameters.gap)[0])
            self.edge_c = parameters.b0 * math.tan(parameters.hgap_m * parameters.fint / parameters.radius_m) / span

    def faces(self):
        result = []
        for row, origin, rotation in self.poses:
            if row.kind == "bend":
                result.append(Face(row.name + " body entrance", origin, rotation[:, 2], False))
                exit_origin = origin + rotation @ np.array([self.p.radius_m * (self.cosine - 1), 0,
                                                             self.p.radius_m * self.sine])
                result.append(Face(row.name + " body exit", exit_origin,
                                   rotation @ np.array([-self.sine, 0, self.cosine]), False))
            elif row.kind == "quadrupole":
                result.append(Face("QACH entrance", origin, rotation[:, 2], True))
                result.append(Face("QACH exit", origin + rotation[:, 2] * self.p.quadrupole_length_m,
                                   rotation[:, 2], False))
        return result

    def __call__(self, position, quad_on=False):
        if abs(position[1]) > 1e-12:
            raise ValueError("The analytic Jacobian is specialized to the y=0 reference")
        field, jacobian = np.zeros(3), np.zeros((3, 3))
        p = self.p
        for row, origin, rotation in self.poses:
            local = rotation.T @ (position - origin)
            if row.kind == "drift":
                continue
            if row.kind == "quadrupole":
                if quad_on:
                    field += rotation @ np.array([0, p.quad_gradient * local[0], 0])
                    jacobian += rotation @ np.array([[0, p.quad_gradient, 0],
                                                     [p.quad_gradient, 0, 0], [0, 0, 0]]) @ rotation.T
                continue
            # GeometryHelper's sector chart and straight entrance/exit extensions.
            exit_dx = local[0] - p.radius_m * (self.cosine - 1)
            exit_dz = local[2] - p.radius_m * self.sine
            exit_z = -self.sine * exit_dx + self.cosine * exit_dz
            if local[2] <= 0:
                z, dz, phi = local[2], np.array([0.0, 0, 1.0]), 0.0
            elif exit_z >= 0:
                z = p.body_length + exit_z
                dz, phi = np.array([-self.sine, 0, self.cosine]), self.theta
            else:
                radius2 = (local[0] + p.radius_m)**2 + local[2]**2
                phi = math.atan2(local[2], local[0] + p.radius_m)
                z = p.radius_m * phi
                dz = p.radius_m * np.array([-local[2], 0, local[0] + p.radius_m]) / radius2
            if z < -p.half_width or z >= p.body_length + p.half_width:
                continue
            entrance, exit_ = enge_first(-z, p.gap), enge_first(z - p.body_length, p.gap)
            if exit_[0] < entrance[0]:
                amplitude, derivative, sign = exit_[0], exit_[1], -1
            else:
                amplitude, derivative, sign = entrance[0], -entrance[1], 1
            field += rotation @ np.array([0, p.b0 * amplitude, 0])
            local_jacobian = np.zeros((3, 3))
            local_jacobian[1, :] = p.b0 * derivative * dz
            local_jacobian[:, 1] = rotation_y(-phi) @ np.array([sign * self.edge_c * derivative,
                                                              0, p.b0 * derivative])
            jacobian += rotation @ local_jacobian @ rotation.T
        return field, jacobian


class VariationalDynamics:
    def __init__(self, field, u0, force_scale):
        self.field = field
        self.u0 = u0
        self.force_scale = force_scale
        self.beta0 = u0 / math.sqrt(1 + u0*u0)

    def values(self, state, quad_on=False):
        k = state[3:6]
        magnitude = np.linalg.norm(k)
        n = k / magnitude
        dn = (np.eye(3) - np.outer(n, n)) / magnitude
        magnetic, gradient = self.field(state[:3], quad_on)
        f = np.zeros(7)
        f[:3] = n
        f[3:6] = self.force_scale * np.cross(n, magnetic)
        small = 1 / self.u0**2
        f[6] = self.beta0 * math.sqrt(magnitude*magnitude + small) / magnitude
        jac = np.zeros((7, 7))
        jac[:3, 3:6] = dn
        jac[3:6, :3] = self.force_scale * cross_matrix(n) @ gradient
        jac[3:6, 3:6] = -self.force_scale * cross_matrix(magnetic) @ dn
        jac[6, 3:6] = -self.beta0 * small * n / (magnitude**2 * math.sqrt(magnitude**2 + small))
        return f, jac

    def rhs(self, quad_on):
        def evaluate(s, augmented):
            f, jac = self.values(augmented[:7], quad_on)
            tangent = augmented[7:].reshape((7, 6))
            return np.concatenate((f, (jac @ tangent).ravel()))
        return evaluate


def initial_state():
    state = np.array([0., 0., 0., 0., 0., 1., 0.])
    tangent = np.zeros((7, 6))
    tangent[0, 0] = tangent[3, 1] = tangent[1, 2] = tangent[4, 3] = tangent[5, 5] = 1
    tangent[6, 4] = -1
    return np.concatenate((state, tangent.ravel()))


def saltation(f_before, f_after, normal):
    surface = np.zeros(7)
    surface[:3] = normal
    denominator = surface @ f_before
    if denominator <= 0:
        raise ValueError("Expected a forward, transverse boundary crossing")
    return np.eye(7) + np.outer(f_after - f_before, surface) / denominator


def encode_exit(dynamics, augmented, quad_on=False):
    state = augmented[:7]
    momentum = np.linalg.norm(state[3:6])
    n = state[3:6] / momentum
    ex = np.array([n[2], 0, -n[0]])  # Bishop frame for this planar trajectory.
    ey = np.array([0., 1., 0.])
    f, _ = dynamics.values(state, quad_on)
    tangent = augmented[7:].reshape((7, 6))
    # Every ray is observed at its own crossing of the same reference-normal plane.
    tangent = tangent - np.outer(f, n @ tangent[:3]) / (n @ f[:3])
    output = np.zeros((6, 7))
    output[0, :3], output[1, 3:6] = ex, ex / momentum
    output[2, :3], output[3, 3:6] = ey, ey / momentum
    beta = dynamics.u0 * momentum / math.sqrt(1 + (dynamics.u0 * momentum)**2)
    output[4, 6] = -beta / dynamics.beta0
    output[5, 3:6] = n / momentum
    return output @ tangent


def integrate(dynamics, stop, faces=(), method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.02):
    augmented = initial_state()
    s = 0.0
    quad_on = False
    records = []
    function_evaluations = 0
    max_energy_error = 0.0
    for face in [*faces, None]:
        event = None
        if face is not None:
            def event(t, z, face=face):
                return (z[:3] - face.origin) @ face.normal
            event.terminal = True
            event.direction = 1
            if event(s, augmented) >= -1e-13:
                raise ValueError(f"Reference is not upstream of {face.name}")
        result = solve_ivp(dynamics.rhs(quad_on), (s, stop), augmented, method=method,
                           rtol=rtol, atol=atol, max_step=max_step, events=event)
        if not result.success:
            raise RuntimeError(result.message)
        function_evaluations += result.nfev
        max_energy_error = max(max_energy_error, float(abs(np.linalg.norm(result.y[3:6], axis=0) - 1).max()))
        s, augmented = float(result.t[-1]), result.y[:, -1]
        if face is not None:
            if result.status != 1:
                raise ValueError(f"Reference did not reach {face.name}")
            f_before, _ = dynamics.values(augmented[:7], quad_on)
            f_after, _ = dynamics.values(augmented[:7], face.next_quad_on)
            tangent = saltation(f_before, f_after, face.normal) @ augmented[7:].reshape((7, 6))
            augmented[7:] = tangent.ravel()
            quad_on = face.next_quad_on
            records.append(dict(face=face.name, s_m=s, position_m=augmented[:3].tolist()))
    matrix = encode_exit(dynamics, augmented, quad_on)
    return matrix, dict(final_state=augmented[:7].tolist(), events=records,
                        function_evaluations=function_evaluations, max_relative_momentum_drift=max_energy_error)


def make_reference(p, **settings):
    u0, mass = particle_constants()
    field = PlanarField(p)
    dynamics = VariationalDynamics(field, u0, p.charge_e * C_LIGHT / (mass * 1e9 * u0))
    return integrate(dynamics, float(layout(p).support_end_m.max()), field.faces(), **settings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "variational-reference-validated")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Use a new output directory to preserve the previous reference")
    args.output.mkdir(parents=True)
    p = load_parameters()
    settings = [dict(name="dop853-coarse", method="DOP853", rtol=1e-10, atol=1e-12, max_step=0.04),
                dict(name="dop853-fine", method="DOP853", rtol=2e-13, atol=2e-15, max_step=0.01),
                dict(name="dop853-finer", method="DOP853", rtol=3e-14, atol=3e-16, max_step=0.005),
                dict(name="rk45-crosscheck", method="RK45", rtol=3e-14, atol=3e-16, max_step=0.005)]
    records, matrices = [], []
    for entry in settings:
        name = entry["name"]
        start = time.monotonic()
        matrix, information = make_reference(p, **{k: v for k, v in entry.items() if k != "name"})
        record = entry | information | dict(elapsed_s=time.monotonic() - start)
        matrices.append(matrix)
        records.append(record)
        np.savetxt(args.output / f"{name}-matrix.txt", matrix, fmt="%.16e")
        print(name, f'{record["elapsed_s"]:.2f} s', "R16/R26", matrix[0, 5], matrix[1, 5], flush=True)
    reference = matrices[2]
    for record, matrix in zip(records, matrices):
        record["max_entry_difference_from_finer"] = float(abs(matrix - reference).max())
    np.savetxt(args.output / "reference-matrix.txt", reference, fmt="%.16e",
               header="Numerically evaluated analytic-field variational reference; not closed form")
    u0, mass = particle_constants()
    j = np.kron(np.eye(3), np.array([[0., 1.], [-1., 0.]]))
    sources = [Path(__file__), ROOT / "analytic_fields.py", ROOT / "validate_opalx_fields.py",
               ROOT / "parameters.json", ROOT / "test_variational_reference.py",
               ROOT.parent / "reference-particle.txt", REPO / "src/Physics/Physics.h"]
    provenance = dict(type="independent analytic-field variational ODE; no OPALX calls or ray differencing",
                      scipy_version=scipy.__version__, u0=u0, electron_mass_GeV=mass,
                      particle_momentum_GeV_c=u0 * mass, field_normalization_momentum_GeV_c=p.momentum_GeV_c,
                      refinements=records, selected_reference="dop853-finer",
                      determinant_error=float(abs(np.linalg.det(reference) - 1)),
                      canonical_J_error=float(abs(reference.T @ j @ reference - j).max()),
                      sources_sha256={str(path.relative_to(REPO)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources})
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Successive refinement differences:", [float(abs(a-b).max()) for a, b in zip(matrices[:-1], matrices[1:])])
    print("Reference matrix:\n", reference)


if __name__ == "__main__":
    main()
