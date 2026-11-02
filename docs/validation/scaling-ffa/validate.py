#!/usr/bin/env python3
"""Reproducible scaling FFA ring tracking validation (no legacy executable needed).

The overlapping tanh sectors telescope to a uniform median-plane field. This
provides an independent circular-orbit oracle for actual production tracking.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MASS_GEV = 0.93827208816
C = 299792458.0
RADIUS, FIELD, SECTORS, PARTICLES = 4.42, 0.1, 16, 8


def prepare(path, dt, sign=1, turns=3, placement='pose'):
    path.mkdir(parents=True, exist_ok=True)
    radius, field, sectors, count = RADIUS, FIELD, SECTORS, PARTICLES
    momentum = C * 1e-9 * field * radius
    bg = momentum / MASS_GEV
    length = 2 * math.pi * radius / sectors
    # A common longitudinal shift separates the injection from the return plane.
    rows = [f'{(i - 3.5) * 1e-5:.17g} 0 0 0 0 {bg:.17g}' for i in range(count)]
    (path / 'particles.dat').write_text(f'{count}\nx px y py z pz\n' + '\n'.join(rows) + '\n')
    lines = ['OPTION, ECHO=FALSE, ENABLEHDF5=FALSE, PSDUMPFREQ=0, STATDUMPFREQ=10, '
             'ENABLELINEARTRANSFERMAPS=FALSE, REPARTFREQ=0, VERSION=10900;']
    for i in range(sectors):
        angle = 2 * math.pi * i / sectors
        pose = (f'X={sign * radius * (math.cos(angle) - 1):.17g}, Y=0, '
                f'Z={radius * math.sin(angle):.17g}, THETA={-sign * angle:.17g}, PHI=0, PSI=0'
                if placement == 'pose' else f'ELEMEDGE={radius * angle:.17g}')
        lines.append(f'M{i}: SCALINGFFAMAGNET, L={length:.17g}, '
                     f'R0={sign * radius}, B0={sign * field}, FIELD_INDEX=0, TAN_DELTA=0, '
                     f'MAX_Y_POWER=4, FRINGE_LENGTH=0.02, '
                     f'AZIMUTHAL_EXTENT={length / 2 + 0.4:.17g}, '
                     f'RADIAL_NEG_EXTENT=0.8, RADIAL_POS_EXTENT=1.2, HEIGHT=1, {pose};')
    lines += ['LATTICE: RING=(' + ','.join(f'M{i}' for i in range(sectors)) + ');',
              f'DIST: DISTRIBUTION, TYPE=FROMFILE, FNAME="particles.dat", NPARTDIST={count};',
              'SOURCE: EMISSIONSOURCE, DISTRIBUTION=DIST;',
              'SOURCES: EMISSIONSOURCELIST=(SOURCE);',
              f'BEAM1: BEAM, PARTICLE=PROTON, NALLOC={count}, SOURCES=SOURCES, '
              'BCHARGE=1e-15, CHARGE=1;',
              'FS: FIELDSOLVER, TYPE=NONE, NX=8, NY=8, NZ=8;',
              f'TRACK, LINE=LATTICE, BEAM=BEAM1, DT={{{dt:.17g}}}, MAXSTEPS=100000, ZSTOP=1000;',
              f'RUN, METHOD="PARALLEL", FIELDSOLVER=FS, TURNS={turns};',
              'ENDTRACK;', 'STOP;']
    (path / 'ring.in').write_text('\n'.join(lines) + '\n')
    (path / 'probes.cfg').write_text('OUTPUT probes.csv\nIDS ' + ' '.join(map(str, range(count)))
                                    + '\nPLANE 0 0 0 0 0 0 1\n')
    return radius, bg


def run(executable, path, ranks, execute=True):
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_PROC_BIND='false',
               OPALX_TEST_PASSIVE_RING_PROBES='probes.cfg')
    command = ['mpiexec', '--host', 'localhost:4', '--map-by', 'slot', '--bind-to', 'none',
               '-n', str(ranks), str(executable), '--info', '1', 'ring.in']
    if execute:
        with (path / 'run.log').open('w') as output:
            subprocess.run(command, cwd=path, env=env, stdout=output, stderr=subprocess.STDOUT,
                           check=True, timeout=240)
    assert 'completed directed turn 3' in (path / 'run.log').read_text()
    data = pd.read_csv(path / 'probes.csv')
    crossing = data[data.status == 'crossing'].copy()
    if crossing.empty:
        raise AssertionError('No passive plane crossings')
    return crossing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, default=ROOT / 'build-omp/src/opalx')
    parser.add_argument('--output', type=Path, default=ROOT / 'output/scaling-ffa')
    parser.add_argument('--check-only', action='store_true', help='Recheck saved runs without rerunning OPALX')
    args = parser.parse_args()
    reports, samples = [], {}
    for sign in (1, -1):
        for placement in ('pose', 'elemedge'):
            for ranks, dt in ((1, 2e-10), (1, 1e-10), (2, 1e-10)):
                name = f'{sign:+d}-{placement}-{ranks}rank-{dt:g}'
                path = args.output.resolve() / name
                if args.check_only:
                    radius, bg = RADIUS, C * 1e-9 * FIELD * RADIUS / MASS_GEV
                else:
                    if path.exists():
                        raise FileExistsError(f'Use a new output directory: {path}')
                    radius, bg = prepare(path, dt, sign, placement=placement)
                data = run(args.executable.resolve(), path, ranks, execute=not args.check_only)
                # IPPL IDs are interleaved by creation rank; FROMFILE rows are
                # contiguous rank blocks. Compare the same physical input rows.
                data['input_row'] = (data.particle_id % ranks) * (8 // ranks) + data.particle_id // ranks
                # The reference controls stop time. A particle at the last return
                # can lie on either side of the final endpoint's roundoff bracket.
                # Compare the two complete preceding turns for every particle.
                data = data[data.turn <= 2].sort_values(['input_row', 'turn']).reset_index(drop=True)
                counts = data.groupby('particle_id').size()
                assert len(counts) == 8 and (counts == 2).all(), counts
                assert not data.duplicated(['input_row', 'turn']).any()
                radial_error = ((data.x_m + sign * radius) ** 2 + data.z_m ** 2) ** .5 - radius
                energy_error = (data.px_mc ** 2 + data.py_mc ** 2 + data.pz_mc ** 2) ** .5 - bg
                # Linear plane interpolation has O(dt^2) position error. The
                # deliberately displaced particles have circular-orbit radial
                # offsets <=35 um. Bound includes their displacement and step error.
                assert radial_error.abs().max() < 5e-5
                assert energy_error.abs().max() < 1e-7
                expected_x = (data.input_row - 3.5) * 1e-5
                return_error = (data.x_m - expected_x).abs().max()
                assert return_error < 3e-6
                samples[(sign, placement, ranks, dt)] = data
                report = dict(case=name, crossings=len(data), min_crossings=int(counts.min()),
                              max_radius_error_m=float(radial_error.abs().max()),
                              max_momentum_error_mc=float(energy_error.abs().max()),
                              max_return_error_m=float(return_error))
                if ranks == 1 and dt == 1e-10:
                    coarse = samples[(sign, placement, 1, 2e-10)]
                    coarse_error = (coarse.x_m - (coarse.input_row - 3.5) * 1e-5).abs().max()
                    assert return_error < 0.4 * coarse_error
                    report['coarse_to_fine_return_error_ratio'] = float(coarse_error / return_error)
                if ranks == 2:
                    serial = samples[(sign, placement, 1, dt)]
                    assert data[['input_row', 'turn']].equals(serial[['input_row', 'turn']])
                    columns = ['time_s', 'x_m', 'y_m', 'z_m', 'px_mc', 'py_mc', 'pz_mc']
                    delta = (data[columns] - serial[columns]).abs().max()
                    assert delta[['x_m', 'y_m', 'z_m', 'px_mc', 'py_mc', 'pz_mc']].max() < 1e-10, delta
                    assert delta.time_s < 1e-16, delta
                    report['max_serial_mpi_position_difference_m'] = float(delta[['x_m', 'y_m', 'z_m']].max())
                    report['max_serial_mpi_momentum_difference_mc'] = float(delta[['px_mc', 'py_mc', 'pz_mc']].max())
                    report['max_serial_mpi_time_difference_s'] = float(delta.time_s)
                if placement == 'elemedge':
                    posed = samples[(sign, 'pose', ranks, dt)]
                    delta = (data[['x_m', 'y_m', 'z_m']] - posed[['x_m', 'y_m', 'z_m']]).abs().max().max()
                    assert delta < 1e-10
                    report['max_placement_position_difference_m'] = float(delta)
                reports.append(report)
                print(json.dumps(report), flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'summary.json').write_text(json.dumps(reports, indent=2) + '\n')


if __name__ == '__main__':
    main()
