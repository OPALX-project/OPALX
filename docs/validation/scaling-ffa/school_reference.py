#!/usr/bin/env python3
"""Independent median-plane polar ODE reference for the school closed orbit.

No OPALX field evaluator or tracking code is used. The 16-fold periodic orbit is
found from one cell. This is a cross-check, not a replacement for the 4D COF.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.special import expit

from school_cof import R0, MASS, ENERGY, CELLS, SPIRAL, K, BF, BD

C = 299792458.
P0 = math.sqrt(ENERGY * (ENERGY + 2 * MASS)) / MASS
PERIOD = 2 * math.pi / CELLS
ENTRY, HALF, STRENGTH = [], [], []
for cell in range(CELLS):
    for start, length, field in [(14.23, 4.5, -BF), (18.73 + 1.52, 2.25, -BD)]:
        ENTRY.append(math.radians(-6 + cell * 22.5 + start - 3 * length / 8))
        HALF.append(math.radians(length) * R0 / 2)
        STRENGTH.append(field)
ENTRY, HALF, STRENGTH = map(np.array, (ENTRY, HALF, STRENGTH))


def field(radius, phi):
    local_phi = (phi - ENTRY + np.pi) % (2 * np.pi) - np.pi
    psi = local_phi - math.tan(SPIRAL) * math.log(radius / R0) - HALF / R0
    u0, u1 = (-R0 * psi - HALF) / .12, (R0 * psi - HALF) / .12
    polynomial = lambda u: -1.3 + u * (3.42 + u * (-1.2 + .23 * u))
    fringe = -1 + expit(-polynomial(u0)) + expit(-polynomial(u1))
    fringe = np.where(np.abs(psi) <= PERIOD, fringe, 0.)
    return (radius / R0) ** K * np.dot(STRENGTH, fringe)


def trace(initial, rtol=5e-12, atol=5e-14, cells=1):
    def rhs(phi, state):
        radius, radial = state
        tangent = math.sqrt(P0 ** 2 - radial ** 2)
        return [radius * radial / tangent,
                tangent - C / (MASS * 1e9) * radius * field(radius, phi)]
    result = solve_ivp(rhs, (0, PERIOD * cells), initial, method='DOP853',
                       rtol=rtol, atol=atol, max_step=PERIOD / 100, dense_output=True)
    if not result.success:
        raise RuntimeError(result.message)
    return result


def solve():
    def residual(initial):
        return (trace(initial).y[:, -1] - initial) / np.array([1, P0])
    result = root(residual, [4.395245598494839, .016268421156887713], tol=1e-10)
    orbit = trace(result.x)
    residual = orbit.y[:, -1] - result.x
    if max(np.abs(residual)) >= 1e-11:
        raise AssertionError((result.message, residual))
    full = trace(result.x, cells=CELLS)
    answer = dict(method='independent SciPy DOP853 polar ODE, one-cell fixed point',
                  radius_m=float(result.x[0]), radial_momentum_mc=float(result.x[1]),
                  x_m=float(result.x[0] - R0), y_m=0, py_mc=0,
                  cell_residual=residual.tolist(),
                  full_turn_residual=(full.y[:, -1] - result.x).tolist(),
                  rtol=5e-12, atol=5e-14)
    return answer, orbit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--energy-mev', type=float, default=12.)
    args = parser.parse_args()
    global P0
    energy = args.energy_mev * 1e-3
    P0 = math.sqrt(energy * (energy + 2 * MASS)) / MASS
    if args.output.exists():
        raise FileExistsError(args.output)
    answer, orbit = solve()
    answer['energy_MeV'] = args.energy_mev
    args.output.mkdir(parents=True)
    (args.output / 'reference.json').write_text(json.dumps(answer, indent=2) + '\n')
    phi = np.linspace(0, PERIOD, 401)
    radius, radial = orbit.sol(phi)
    np.savetxt(args.output / 'reference-cell.csv', np.column_stack([phi, radius, radial]),
               delimiter=',', header='phi_rad,radius_m,radial_momentum_mc', comments='')
    print(json.dumps(answer, indent=2))


if __name__ == '__main__':
    main()
