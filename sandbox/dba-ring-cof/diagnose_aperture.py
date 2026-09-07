#!/usr/bin/env python3
"""Geometry-only reproducer of the COF all-element aperture false positive.

No tracking, aperture changes, or production-code changes. Uses the same
explicit entrance poses and Y-axis rotation convention as the input deck.
ElementBase::applyToReferenceParticle tests 0 <= local z < L followed by
transverse aperture exclusion, without ring-path ownership information.
"""
from pathlib import Path
import runpy
import numpy as np

ROOT = Path(__file__).resolve().parent


def main():
    model = runpy.run_path(str(ROOT.parent / 'map-3/ring/build_ring.py'))['geometry'](0)
    elements = {row['name']: row for row in model['elements']}
    drift = elements['C1_D1']
    theta = drift['theta_rad']
    point = np.array([drift['x_m'] + .5*np.sin(theta), 0.,
                      drift['z_m'] + .5*np.cos(theta)])
    print('Exact design orbit, midpoint C1_D1; global XYZ [m]:', point)
    for name in ('C1_D1', 'C4_D2'):
        row = elements[name]
        theta = row['theta_rad']
        delta = point - np.array([row['x_m'], 0., row['z_m']])
        local = np.array([np.cos(theta)*delta[0]-np.sin(theta)*delta[2],
                          delta[1], np.sin(theta)*delta[0]+np.cos(theta)*delta[2]])
        in_z = 0 <= local[2] < row['length_m']
        print(name, 'local XYZ [m]:', local, 'longitudinal test:', in_z)
        assert in_z
        if name == 'C1_D1':
            assert abs(local[0]) < 1e-14
        else:
            assert abs(local[0]) > 7.8
    print('Confirmed: an on-orbit point is in the opposite drift longitudinal slab,')
    print('but far outside its aperture. Checking all elements causes a false loss.')


if __name__ == '__main__':
    main()
