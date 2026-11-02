#!/usr/bin/env python3
"""Plot the independently checked COF orbit and a production reference traversal."""
import argparse
from io import StringIO
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from school_cof import R0
from school_reference import solve, PERIOD


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--track-dir', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(args.output.resolve() / '.matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(args.output.resolve() / '.cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.labelsize': 10, 'axes.titlesize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'savefig.dpi': 220})
    solution = json.loads((args.track_dir / 'orbit.json').read_text())
    text = (args.track_dir / 'data/school_DesignPath.dat').read_text()
    numeric = '\n'.join(line.split('\t', 1)[0] for line in text.splitlines())
    path = pd.read_csv(StringIO(numeric), sep=r'\s+', comment='#', header=None)
    theta = np.unwrap(np.arctan2(path[3].to_numpy(), path[1].to_numpy() + R0))
    radii = np.hypot(path[1].to_numpy() + R0, path[3].to_numpy())
    reference, cell = solve()
    phi = np.linspace(0, PERIOD, 401)
    r = cell.sol(phi)[0]
    full_phi = np.linspace(0, 2 * np.pi, 6401)
    full_r = cell.sol(full_phi % PERIOD)[0]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3), layout='constrained')
    ax = axes[0]
    ax.plot(R0 * np.cos(full_phi), R0 * np.sin(full_phi), color='.65', ls='--', lw=1,
            label='Nominal radius: 4.42 m')
    ax.plot(full_r * np.cos(full_phi), full_r * np.sin(full_phi), color='#165b9b', lw=1.5,
            label='Independent closed orbit')
    ax.plot(path[1] + R0, path[3], color='#b95026', lw=.8, alpha=.8,
            label='Production Boris reference')
    ax.scatter([R0 + solution['coordinates'][0]], [0], s=27, color='black', zorder=4)
    ax.set(xlabel='Horizontal X from ring centre [m]', ylabel='Horizontal Z [m]',
           title='16-cell school FFA · 12 MeV protons', aspect='equal')
    ax.legend(loc='center', fontsize=8)
    ax = axes[1]
    ax.axvspan(6.5425, 11.0425, color='#4c8bb8', alpha=.16, label='F profile span')
    ax.axvspan(13.40625, 15.65625, color='#df8951', alpha=.22, label='D profile span')
    ax.plot(np.degrees(phi), r, color='#165b9b', lw=1.8, label='Independent polar ODE')
    selected = (theta >= -1e-12) & (theta <= PERIOD)
    ax.plot(np.degrees(theta[selected]), radii[selected], color='#b95026', ls='--', lw=1,
            label='Production Boris reference')
    ax.set(xlabel='Azimuth from the launch section [deg]', ylabel='Radius [m]',
           xlim=(0, 22.5), title=f'One periodic cell · launch radius {R0 + solution["coordinates"][0]:.9f} m')
    ax.grid(axis='y', alpha=.2)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.20), fontsize=8, ncol=2)
    for extension in ['png', 'pdf']:
        fig.savefig(args.output / f'school-closed-orbit.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    main()
