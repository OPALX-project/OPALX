#!/usr/bin/env python3
"""Generate and solve the original 16-cell, 12 MeV FFA25 school lattice.

Legacy F and D placement is sequential: D's start is measured after F_END.
B0 changes sign under the old (x,y horizontal; z vertical) to new (x,z
horizontal; y vertical) reflection. Legacy Enge.getCentreLength()=x0/2
also requires a -3*L/8 entrance shift to preserve the old field centres.
Zero-field circular design arcs fill nominal
geometrical gaps; they add no magnetic field.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
R0, MASS, ENERGY, CELLS = 4.42, 0.938272, 0.012, 16
SPIRAL = 3.142 / 4  # Exact value passed by the original notebook, not pi/4.
K = 7.717022064259446
BF, BD = -0.6532097983351404, 0.502765236428169


def lattice():
    lines = ['// FFA25_OPAL.ipynb, cell 3: full 16-cell ring, SI units.',
             '// Legacy MAGNET_END advances the next element; preserve -6 deg phase.']
    segments = []
    physical = []
    for cell in range(CELLS):
        start = -6 + cell * 22.5
        for kind, entry, length, field in [('F', start + 14.23, 4.5, -BF),
                                           ('D', start + 18.73 + 1.52, 2.25, -BD)]:
            name = f'C{cell:02d}_{kind}'
            # Old Enge centre length was x0/2=L/4; setupEndField added
            # half of that (L/8) to MAGNET_START. Current OPALX adds L/2.
            # Translate the entrance by -3L/8 to retain the actual old field.
            entry -= 3 * length / 8
            physical.append((entry, entry + length, name, field))
            half = math.radians(length) * R0 / 2
            lines.append(f'{name}_END: ASYMMETRIC_ENGE, '
                         f'START_HALF_LENGTH={half:.17g}, END_HALF_LENGTH={half:.17g}, '
                         'START_FRINGE_LENGTH=0.12, END_FRINGE_LENGTH=0.12, '
                         'START_COEFFICIENTS={-1.3,3.42,-1.2,0.23}, '
                         'END_COEFFICIENTS={-1.3,3.42,-1.2,0.23};')
    cursor = 0.0
    for entry, end, name, field in sorted(physical):
        if entry > cursor:
            segments.append((cursor, entry, f'GAP{len(segments):02d}', 0, False))
        segments.append((entry, end, name, field, True))
        cursor = end
    segments.append((cursor, 360., 'GAP_FINAL', 0, False))
    names = []
    for entry, end, name, field, real_magnet in segments:
        phi = math.radians(entry)
        length = math.radians(end - entry) * R0
        model = f'END_FIELD_MODEL={name}_END' if real_magnet else f'L={length:.17g}, FRINGE_LENGTH=0.12'
        extent = math.radians(22.5) * R0 if real_magnet else length / 2 + 0.12
        lines.append(f'{name}: SCALINGFFAMAGNET, {model}, R0={R0}, B0={field:.17g}, '
                     f'FIELD_INDEX={K}, TAN_DELTA={math.tan(SPIRAL):.17g}, MAX_Y_POWER=4, '
                     f'RADIAL_NEG_EXTENT=0.8, RADIAL_POS_EXTENT=1.2, HEIGHT=1, '
                     f'AZIMUTHAL_EXTENT={extent:.17g}, '
                     f'X={R0 * (math.cos(phi) - 1):.17g}, Y=0, '
                     f'Z={R0 * math.sin(phi):.17g}, THETA={-phi:.17g}, PHI=0, PSI=0;')
        names.append(name)
    lines.append('SCHOOL: RING=(' + ','.join(names) + ');')
    return lines


def prepare(path, dt=1e-10, fd=1e-5, initial=(-0.024754401505161, 0.016268421156888, 0., 0.), track_turns=0,
            integrator='DOP853', track_dt=1e-10, energy_mev=12.):
    path.mkdir(parents=True, exist_ok=True)
    energy = energy_mev * 1e-3
    p0 = math.sqrt(energy * (energy + 2 * MASS))
    bg = p0 / MASS
    lines = ['OPTION, ECHO=FALSE, ENABLEHDF5=FALSE, ENABLELINEARTRANSFERMAPS=FALSE, '
             'PSDUMPFREQ=0, STATDUMPFREQ=100, REPARTFREQ=0, VERSION=10900;']
    lines += lattice()
    lines += ['DIST: DISTRIBUTION, TYPE=GAUSS, NPARTDIST=8, '
              'SIGMAX=0, SIGMAY=0, SIGMAZ=0, SIGMAPX=0, SIGMAPY=0, SIGMAPZ=0;',
              'SOURCE: EMISSIONSOURCE, DISTRIBUTION=DIST;',
              'SOURCES: EMISSIONSOURCELIST=(SOURCE);',
              f'BEAM1: BEAM, PARTICLE=PROTON, MASS={MASS}, PC={p0:.17g}, '
              'NALLOC=8, SOURCES=SOURCES, BCHARGE=1e-15, CHARGE=1;',
              f'ORBIT: COF, LINE=SCHOOL, BEAM=BEAM1, DT={dt:.17g}, '
              f'TIMEINTEGRATOR="{integrator}", MAXPATH=40, MAXSTEPS=100000, '
              f'X={initial[0]:.17g}, PX={initial[1]:.17g}, Y={initial[2]:.17g}, PY={initial[3]:.17g}, '
              f'FDSTEP={{{fd:.17g},{fd * bg:.17g},{fd:.17g},{fd * bg:.17g}}}, '
              f'SCALES={{1,{bg:.17g},1,{bg:.17g}}}, '
              'XTOL=1e-10, PTOL=1e-11, MAXIT=20, OUTPUT="orbit.json";']
    if track_turns:
        lines += ['FS: FIELDSOLVER, TYPE=NONE, NX=8, NY=8, NZ=8;',
                  f'TRACK, LINE=SCHOOL, BEAM=BEAM1, INITIALORBIT=ORBIT, DT={{{track_dt:.17g}}}, '
                  'MAXSTEPS=500000, ZSTOP=1000;',
                  f'RUN, METHOD="PARALLEL", FIELDSOLVER=FS, TURNS={track_turns};', 'ENDTRACK;']
        (path / 'probes.cfg').write_text('OUTPUT probes.csv\nIDS 0 1 2 3 4 5 6 7\n'
                                        'PLANE 0 0 0 0 0 0 1\n')
    lines.append('STOP;')
    (path / 'school.in').write_text('\n'.join(lines) + '\n')


def run(path, executable, ranks=1, track_turns=0, energy_mev=12.):
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_PROC_BIND='false')
    if track_turns:
        env['OPALX_TEST_PASSIVE_RING_PROBES'] = 'probes.cfg'
    else:
        env.pop('OPALX_TEST_PASSIVE_RING_PROBES', None)
    with (path / 'run.log').open('w') as stream:
        subprocess.run(['mpiexec', '--host', 'localhost:4', '--map-by', 'slot', '--bind-to', 'none',
                        '-n', str(ranks), str(executable), '--info', '1', 'school.in'],
                       cwd=path, env=env, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=600)
    result = json.loads((path / 'orbit.json').read_text())
    assert result['converged'] and abs(result['energy_MeV'] - energy_mev) < 1e-9
    assert max(abs(result['residual'][i]) for i in (0, 2)) < 1e-10
    assert max(abs(result['residual'][i]) for i in (1, 3)) < 1e-11
    if track_turns:
        assert f'completed directed turn {track_turns},' in (path / 'run.log').read_text()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--executable', type=Path, default=ROOT / 'build-omp/src/opalx')
    parser.add_argument('--dt', type=float, default=1e-10)
    parser.add_argument('--fdstep', type=float, default=1e-5)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--track-turns', type=int, default=0)
    parser.add_argument('--track-dt', type=float, default=1e-10)
    parser.add_argument('--integrator', choices=['DOP853', 'RK4', 'BORIS'], default='DOP853')
    parser.add_argument('--energy-mev', type=float, default=12.)
    parser.add_argument('--initial', type=float, nargs=4, default=[-0.024754401505161, 0.016268421156888, 0., 0.])
    args = parser.parse_args()
    path = args.output.resolve()
    if path.exists():
        raise FileExistsError(f'Use a new output directory: {path}')
    prepare(path, args.dt, args.fdstep, args.initial, args.track_turns, args.integrator, args.track_dt, args.energy_mev)
    result = run(path, args.executable.resolve(), args.ranks, args.track_turns, args.energy_mev)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
