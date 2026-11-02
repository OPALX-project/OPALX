#!/usr/bin/env python3
"""School FFA tunes from production particle tracking at equivalent cell planes.

Use a tracked closed reference and deterministic fixed-energy +/- displacements.
Complex Courant-Snyder signals distinguish f from 1-f; the integer ring tune is
resolved through sixteen cell samples, rather than a once-per-turn FFT.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from school_cof import lattice, R0, MASS, CELLS, ROOT
from school_linear_optics import optics


def peak(signal):
    """Hann-windowed complex spectrum, followed by local continuous refinement."""
    z = np.asarray(signal, dtype=complex)
    if len(z) < 64 or not np.isfinite(z).all() or np.std(z) == 0:
        raise ValueError('Need at least 64 finite, nonconstant complex samples')
    z = (z - np.mean(z)) * np.hanning(len(z))
    count = 16 * len(z)
    power = abs(np.fft.fft(z, n=count)) ** 2
    frequencies = np.arange(count) / count
    # Discard only the DC neighbourhood, retaining both frequency orientations.
    power[(frequencies < 2 / len(z)) | (frequencies > 1 - 2 / len(z))] = 0
    i = int(power.argmax())
    centre = frequencies[i]
    indices = np.arange(len(z))
    def objective(frequency):
        return -abs(np.dot(z, np.exp(-2j * np.pi * frequency * indices))) ** 2
    fit = minimize_scalar(objective, bounds=(centre - 1 / count, centre + 1 / count),
                          method='bounded', options={'xatol': 1e-14})
    if not fit.success:
        raise RuntimeError(fit.message)
    return float(fit.x % 1)


def prepare(path, dt, turns, amplitude, ranks):
    path.mkdir(parents=True, exist_ok=False)
    orbit = json.loads((ROOT / 'docs/validation/scaling-ffa/school-closed-orbit.json').read_text())
    position = np.array(orbit['position_m'])
    momentum = np.array(orbit['momentum_mc'])
    offsets = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (0, 0)]
    rows = []
    for radial, vertical in offsets:
        point = position + amplitude * np.array([radial, vertical, 0.])
        rows.append(' '.join(f'{v:.17g}' for v in
                             [point[0], momentum[0], point[1], momentum[1], point[2], momentum[2]]))
    (path / 'particles.dat').write_text('8\nx px y py z pz\n' + '\n'.join(rows) + '\n')
    statements = ['OPTION, ECHO=FALSE, ENABLEHDF5=FALSE, ENABLELINEARTRANSFERMAPS=FALSE, '
                  'PSDUMPFREQ=0, STATDUMPFREQ=1000, REPARTFREQ=0, VERSION=10900;']
    statements += lattice()
    statements += ['DIST: DISTRIBUTION, TYPE=FROMFILE, FNAME="particles.dat", NPARTDIST=8;',
                   'SOURCE: EMISSIONSOURCE, DISTRIBUTION=DIST;',
                   'SOURCES: EMISSIONSOURCELIST=(SOURCE);',
                   f'BEAM1: BEAM, PARTICLE=PROTON, MASS={MASS}, NALLOC=8, SOURCES=SOURCES, '
                   'BCHARGE=1e-15, CHARGE=1;',
                   'FS: FIELDSOLVER, TYPE=NONE, NX=8, NY=8, NZ=8;',
                   f'TRACK, LINE=SCHOOL, BEAM=BEAM1, DT={{{dt:.17g}}}, MAXSTEPS=2000000, ZSTOP=10000;',
                   f'RUN, METHOD="PARALLEL", FIELDSOLVER=FS, TURNS={turns};', 'ENDTRACK;', 'STOP;']
    (path / 'tunes.in').write_text('\n'.join(statements) + '\n')
    probes = ['OUTPUT probes.csv', 'IDS 0 1 2 3 4 5 6 7']
    for plane in range(CELLS):
        theta = 2 * np.pi * plane / CELLS
        probes.append(f'PLANE {plane} {-R0:.17g} 0 0 {-np.sin(theta):.17g} 0 {np.cos(theta):.17g}')
    (path / 'probes.cfg').write_text('\n'.join(probes) + '\n')
    config = dict(dt_s=dt, turns=turns, amplitude_m=amplitude, ranks=ranks, energy_MeV=12.,
                  offsets=offsets, initial_position_m=position.tolist(), initial_momentum_mc=momentum.tolist())
    (path / 'config.json').write_text(json.dumps(config, indent=2) + '\n')


def run(path, executable, ranks):
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_PROC_BIND='false',
               OPALX_TEST_PASSIVE_RING_PROBES='probes.cfg')
    with (path / 'run.log').open('w') as stream:
        subprocess.run(['mpiexec', '--host', 'localhost:4', '--map-by', 'slot', '--bind-to', 'none',
                        '-n', str(ranks), str(executable), '--info', '1', 'tunes.in'],
                       cwd=path, env=env, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=max(900, 30 * json.loads((path / "config.json").read_text())["turns"]), check=True)


def analyze(path, linear=None, window_turns=None):
    config = json.loads((path / 'config.json').read_text())
    if f'completed directed turn {config["turns"]},' not in (path / 'run.log').read_text():
        raise AssertionError('Requested production tracking did not finish')
    analysis_turns = config['turns'] - 2 if window_turns is None else window_turns
    if not 4 <= analysis_turns <= config['turns'] - 2:
        raise ValueError('Analysis window must leave the first/last partial turns out')
    data = pd.read_csv(path / 'probes.csv')
    if not (data.status == 'crossing').all():
        raise AssertionError('A selected particle went missing or a diagnostic failed')
    ranks = config['ranks']
    # FROMFILE distributes contiguous blocks, whereas IPPL IDs interleave ranks.
    data['row'] = (data.particle_id % ranks) * (8 // ranks) + data.particle_id // ranks
    data['cell'] = (CELLS * (data.turn - 1) + data.plane_id
                    + np.where(data.plane_id == 0, CELLS, 0))
    data = data[(data.cell >= CELLS) & (data.cell < CELLS * (analysis_turns + 1))].copy()
    if data.duplicated(['row', 'cell']).any():
        raise AssertionError('Duplicate particle/cell crossing')
    theta = 2 * np.pi * data.plane_id / CELLS
    data['radius'] = np.hypot(data.x_m + R0, data.z_m)
    data['pr'] = data.px_mc * np.cos(theta) + data.pz_mc * np.sin(theta)
    data['pnorm'] = np.sqrt(data.px_mc ** 2 + data.py_mc ** 2 + data.pz_mc ** 2)
    expected = np.arange(CELLS, CELLS * (analysis_turns + 1))
    frames = {}
    for row in range(8):
        frame = data[data.row == row].set_index('cell').sort_index()
        if not np.array_equal(frame.index.to_numpy(), expected):
            raise AssertionError(f'Missing/extra cell samples for input row {row}')
        frames[row] = frame
    reference = frames[0]
    linear = linear or optics()
    result = dict(config=config, samples=len(expected), analysis_turns=analysis_turns,
                  method='complex Courant-Snyder spectrum of cell-plane samples', tunes={})
    for plane, rows, qname, pname in [('horizontal', [1, 2], 'radius', 'pr'),
                                      ('vertical', [3, 4], 'y_m', 'py_mc')]:
        alpha = linear[plane]['alpha']
        beta = linear[plane]['beta_m_per_mc']
        values = []
        for row in rows:
            q = (frames[row][qname] - reference[qname]).to_numpy()
            p = (frames[row][pname] - reference[pname]).to_numpy()
            signal = q - 1j * (alpha * q + beta * p)
            frequency = peak(signal)
            values.append(dict(row=row, cell_tune=frequency, ring_tune=CELLS * frequency,
                               rms_displacement_m=float(np.sqrt(np.mean(q * q)))))
        mean = float(np.mean([v['ring_tune'] for v in values]))
        result['tunes'][plane] = dict(ring_tune=mean, plus_minus_difference=abs(values[0]['ring_tune'] - values[1]['ring_tune']),
                                     linear_ring_tune=linear[plane]['ring_tune'],
                                     linear_difference=mean - linear[plane]['ring_tune'], particles=values)
    p0 = np.linalg.norm(config['initial_momentum_mc'])
    result['max_relative_momentum_norm_drift'] = float(abs(data.pnorm - p0).max() / p0)
    duplicate = frames[7][['radius', 'pr', 'y_m', 'py_mc']] - reference[['radius', 'pr', 'y_m', 'py_mc']]
    result['reference_duplicate_max_difference'] = float(abs(duplicate).max().max())
    (path / ('tunes.json' if window_turns is None else f'tunes-{window_turns}turns.json')).write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--executable', type=Path, default=ROOT / 'build-omp/src/opalx')
    parser.add_argument('--turns', type=int, default=50)
    parser.add_argument('--dt', type=float, default=1e-10)
    parser.add_argument('--amplitude', type=float, default=1e-4)
    parser.add_argument('--ranks', type=int, choices=[1, 2], default=1)
    parser.add_argument('--window-turns', type=int)
    parser.add_argument('--analyze-only', action='store_true')
    args = parser.parse_args()
    path = args.output.resolve()
    if not args.analyze_only:
        if args.turns < 8 or args.dt <= 0 or args.amplitude <= 0:
            raise ValueError('Need >=8 turns and positive timestep/amplitude')
        prepare(path, args.dt, args.turns, args.amplitude, args.ranks)
        run(path, args.executable.resolve(), args.ranks)
    print(json.dumps(analyze(path, window_turns=args.window_turns), indent=2))


if __name__ == '__main__':
    main()
