#!/usr/bin/env python3
"""Reproducible hard-edge optics scan; no OPALX tracking or field retuning.

Coordinates are (x,x',y,y',zeta,delta), lengths m, K1 m^-2. The
electron K1 convention is inherited from the validated map-2 reference.
The two-family variant replaces 0.2 m of the first drift by a quadrupole,
leaving all bend poses and the circumference unchanged. No achromat constraint.
"""
import json
from pathlib import Path
import runpy

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OPTICS = runpy.run_path(str(ROOT.parent / 'map-2/check_maps.py'))


def cell(k1, k2=None):
    bend = OPTICS['sector_bend'](2., np.pi / 6)
    drift = OPTICS['drift']
    def quad(k, length):
        return drift(length) if k == 0 else OPTICS['quadrupole'](k, length)
    first = [drift(1.)] if k2 is None else [drift(.4), quad(k2, .2), drift(.4)]
    return OPTICS['product'](bend, *first, quad(k1, .2), drift(1.), bend)


def describe(k1, k2=None):
    m = cell(k1, k2)
    h = np.array([np.trace(m[:2, :2]), np.trace(m[2:4, 2:4])]) / 2
    stable = bool(np.all(np.abs(h) < 1))
    row = dict(k1_m2=float(k1), k2_m2=None if k2 is None else float(k2),
               half_trace_x=float(h[0]), half_trace_y=float(h[1]), stable=stable)
    if stable:
        # Principal conjugate-pair phases only, not integer/branch assignments.
        phase = (6 * np.arccos(h) / (2*np.pi)) % 1
        q = np.minimum(phase, 1-phase)
        row.update(qx_principal=float(q[0]), qy_principal=float(q[1]))
    return row


def main():
    output = ROOT / 'analytic-results'
    output.mkdir(exist_ok=True)
    single = pd.DataFrame([describe(k) for k in np.linspace(-15, 15, 6001)])
    single.to_csv(output / 'one-family.csv', index=False)
    candidates = single[single.stable].copy()
    double = None
    if candidates.empty:
        double = pd.DataFrame([describe(k1, k2) for k1 in np.linspace(-10, 10, 161)
                               for k2 in np.linspace(-10, 10, 161)])
        double.to_csv(output / 'two-family.csv', index=False)
        candidates = double[double.stable].copy()
    # Target nonresonant, separated principal phases; retain cell stability margin.
    candidates = candidates[(candidates.half_trace_x.abs() < .95)
                            & (candidates.half_trace_y.abs() < .95)]
    candidates['score'] = ((candidates.qx_principal-.23)**2
                           + (candidates.qy_principal-.31)**2)
    best = candidates.sort_values(['score', 'k1_m2', 'k2_m2']).iloc[0]
    k2 = None if pd.isna(best.k2_m2) else float(best.k2_m2)
    m = cell(best.k1_m2, k2)
    ring = np.linalg.matrix_power(m, 6)
    np.savetxt(output / 'cell-slope-map.txt', m, fmt='%.16e')
    np.savetxt(output / 'ring-slope-map.txt', ring, fmt='%.16e')
    p0 = OPTICS['P']  # p/(m_e c), not GeV/c
    transform = np.diag([1., p0, 1., p0])
    mechanical = transform @ ring[:4, :4] @ np.linalg.inv(transform)
    np.savetxt(output / 'ring-mechanical-map.txt', mechanical, fmt='%.16e')
    dispersion = np.linalg.solve(np.eye(2)-ring[:2, :2], ring[:2, 5])
    result = dict(one_family_stable_samples=int(single.stable.sum()),
                  two_family_stable_samples=None if double is None else int(double.stable.sum()),
                  selected=best.to_dict(), circumference_m=6*(2*np.pi/3+2.2),
                  momentum_GeV_c=.2505104781131461, normalized_momentum=p0,
                  periodic_dispersion_m=float(dispersion[0]),
                  periodic_dispersion_slope=float(dispersion[1]),
                  ring_R16_m=float(ring[0, 5]), ring_R26=float(ring[1, 5]),
                  scope='Analytic hard-edge linear optics only; COF not yet validated')
    (output / 'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
