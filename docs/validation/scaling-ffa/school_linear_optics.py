#!/usr/bin/env python3
"""Independent one-cell linear optics for the fixed-12 MeV school orbit."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

import school_reference as ref
from school_cof import R0, MASS, SPIRAL, K, CELLS


def field_and_gradients(radius, phi):
    local = (phi - ref.ENTRY + np.pi) % (2 * np.pi) - np.pi
    psi = local - np.tan(SPIRAL) * np.log(radius / R0) - ref.HALF / R0
    u0, u1 = (-R0 * psi - ref.HALF) / .12, (R0 * psi - ref.HALF) / .12
    def edge(u):
        e = expit(-(-1.3 + u * (3.42 + u * (-1.2 + .23 * u))))
        derivative = -e * (1 - e) * (3.42 + u * (-2.4 + .69 * u))
        return e, derivative
    e0, d0 = edge(u0)
    e1, d1 = edge(u1)
    mask = np.abs(psi) <= ref.PERIOD
    scale = (radius / R0) ** K
    field = scale * np.dot(ref.STRENGTH, np.where(mask, -1 + e0 + e1, 0.))
    dphi = scale * np.dot(ref.STRENGTH, np.where(mask, (d1 - d0) * R0 / .12, 0.))
    dr = (K * field - np.tan(SPIRAL) * dphi) / radius
    return field, dr, dphi


def optics(rtol=2e-11, max_step=None):
    orbit, _ = ref.solve()
    initial = np.r_[orbit['radius_m'], orbit['radial_momentum_mc'], np.eye(2).ravel(), np.eye(2).ravel()]
    strength = ref.C / (MASS * 1e9)
    def rhs(phi, state):
        radius, radial = state[:2]
        tangent = np.sqrt(ref.P0 ** 2 - radial ** 2)
        field, dr, dphi = field_and_gradients(radius, phi)
        horizontal = np.array([[radial / tangent, radius * ref.P0 ** 2 / tangent ** 3],
                               [-strength * (field + radius * dr), -radial / tangent]])
        vertical = np.array([[0., radius / tangent],
                             [strength * (radius * dr - radial / tangent * dphi), 0.]])
        return np.r_[radius * radial / tangent, tangent - strength * radius * field,
                     (horizontal @ state[2:6].reshape(2, 2)).ravel(),
                     (vertical @ state[6:10].reshape(2, 2)).ravel()]
    result = solve_ivp(rhs, (0, ref.PERIOD), initial, method='DOP853', rtol=rtol,
                       atol=rtol * 1e-2, max_step=max_step or ref.PERIOD / 200)
    if not result.success:
        raise RuntimeError(result.message)
    answer = dict(energy_MeV=12., cells=CELLS, method='independent polar variational equations',
                  rtol=rtol, max_step_rad=max_step or ref.PERIOD / 200)
    for name, block in [('horizontal', slice(2, 6)), ('vertical', slice(6, 10))]:
        history = result.y[block].T.reshape(-1, 2, 2)
        matrix = history[-1]
        cosine = np.trace(matrix) / 2
        if abs(cosine) >= 1:
            raise AssertionError(f'Unstable cell: {name}, trace/2={cosine}')
        sine = np.copysign(np.sqrt(1 - cosine ** 2), matrix[0, 1])
        phase = np.arctan2(sine, cosine) % (2 * np.pi)
        beta = matrix[0, 1] / sine
        if beta <= 0 or abs(np.linalg.det(matrix) - 1) > 1e-9:
            raise AssertionError((name, matrix))
        alpha = (matrix[0, 0] - matrix[1, 1]) / (2 * sine)
        q = history[:, 0, 0] - alpha / beta * history[:, 0, 1]
        p = history[:, 1, 0] - alpha / beta * history[:, 1, 1]
        continuous_phase = np.unwrap(np.angle(q - 1j * (alpha * q + beta * p)))[-1]
        # The trace only determines phase modulo 2pi. Follow a nonzero solution
        # continuously through the cell to establish its integer winding.
        if abs(np.remainder(continuous_phase - phase + np.pi, 2 * np.pi) - np.pi) > 1e-9:
            raise AssertionError((name, continuous_phase, phase))
        phase = continuous_phase
        answer[name] = dict(matrix=matrix.tolist(), cell_phase_rad=float(phase),
                            cell_tune=float(phase / (2 * np.pi)),
                            ring_tune=float(CELLS * phase / (2 * np.pi)),
                            alpha=float(alpha), phase_samples=len(history),
                            beta_m_per_mc=float(beta), determinant=float(np.linalg.det(matrix)))
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = optics()
    refined = optics(rtol=1e-12, max_step=ref.PERIOD / 400)
    for plane in ['horizontal', 'vertical']:
        difference = abs(result[plane]['ring_tune'] - refined[plane]['ring_tune'])
        if difference > 1e-8:
            raise AssertionError((plane, difference))
        result[plane]['refinement_difference'] = difference
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
