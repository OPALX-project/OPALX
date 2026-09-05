#!/usr/bin/env python3
"""Reproduce the first accelerating turn using installed OPAL 2022.1, one rank.

Generate isolated inputs from cyclotron2.in; preserve the supplied original.
The named trim-coil conversion follows the supplied modern cyclotron1.in.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps-per-turn', type=int, default=2880)
    parser.add_argument('--integrator', choices=['LF2', 'RK4'], default='LF2')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.steps_per_turn < 720:
        parser.error('Use at least 720 steps per nominal turn')
    work = args.output.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = (root/'opal/cyclotron2.in').read_text()
    source = re.sub(r'(?m)^(\w+\s*=[^;\n]+;)[ \t]*$', r'REAL \1', source)
    source = source.replace('REAL FREQ=', 'FREQ=')
    source = source.replace('OPTION, PSDUMPLOCALFRAME= TRUE;',
                            'OPTION, VERSION=20200, SPTDUMPFREQ=1, PSDUMPFRAME=GLOBAL;')
    source = source.replace('TCR1=4350, TCR2=4470, MBTC=14e-3, SLPTC=6.0', 'TRIMCOIL=TC1')
    source = source.replace('Ring: CYCLOTRON',
        'TC1: TRIMCOIL, TYPE="PSI-BFIELD-MIRRORED", RMIN=4350, RMAX=4470, '
        'BMAX=1.4e-3, SLPTC=0.6;\nRing: CYCLOTRON')
    source = source.replace('DISTRIBUTION=FROMFILE', 'TYPE=FROMFILE')
    source = source.replace('Ring: CYCLOTRON', 'RingCycl: CYCLOTRON')
    source = source.replace('(Ring,', '(RingCycl,')
    for name in ['bfield.dat', 'rffield1.dat', 'rffield2.dat', 'dist2.dat']:
        source = source.replace(f'"{name}"', f'"{root / "opal" / name}"')
    # Small margin brackets the actual directed turn, not only the nominal RF period.
    source = source.replace('MAXSTEPS= 720*15, STEPSPERTURN= 720',
        f'MAXSTEPS={int(args.steps_per_turn*1.1)}, STEPSPERTURN={args.steps_per_turn}, '
        f'TIMEINTEGRATOR="{args.integrator}"')
    (work/'acceleration.in').write_text(source)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    with (work/'run.log').open('w') as log:
        subprocess.run(['bash', '-c',
            'source /Users/adelmann/OPAL-2022.1/etc/profile.d/opal.sh && '
            '/Users/adelmann/OPAL-2022.1/bin/opal acceleration.in'],
            cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    orbit = pd.read_csv(work/'acceleration-trackOrbit.dat', sep=r'\s+', comment='#',
                        names=['id', 'x', 'px', 'y', 'py', 'z', 'pz'])
    r = orbit[['x', 'y', 'z']].to_numpy()
    p = orbit[['px', 'py', 'pz']].to_numpy()
    normal = p[0]/np.linalg.norm(p[0])
    distance = (r-r[0])@normal
    armed = False
    crossing = None
    for i in range(1, len(r)):
        if distance[i] < -1e-9:
            armed = True
        if armed and distance[i] >= 0 and p[i]@normal > 0:
            crossing = i
            break
    if crossing is None:
        raise RuntimeError('No directed first-turn crossing in the bounded run')
    dt = 6/(50.65e6*args.steps_per_turn)
    # Use the launch energy to infer the exact mass convention from the rounded dump.
    # Energy uncertainty is therefore limited by legacy ASCII momentum precision.
    mass = 72/(np.sqrt(1+np.dot(p[0], p[0]))-1)
    orbit['time_s'] = np.arange(len(orbit))*dt
    orbit['kinetic_MeV'] = (np.sqrt(1+np.sum(p*p, axis=1))-1)*mass
    orbit.iloc[:crossing+1].to_csv(work/'first-turn.csv', index=False)
    changes = np.diff(orbit.kinetic_MeV.to_numpy()[:crossing+1])
    event_steps = np.flatnonzero(np.abs(changes) > 1e-3)+1
    if len(event_steps) != 5:
        raise RuntimeError(f'Expected five first-turn RF kicks, found {len(event_steps)}')
    log = (work/'run.log').read_text()
    cavity_lines = [line for line in log.splitlines() if 'transit time factor=' in line]
    metrics = dict(integrator=args.integrator, ranks=1,
        steps_per_nominal_turn=args.steps_per_turn, dt_s=dt,
        directed_turn_step=crossing, directed_turn_time_s=crossing*dt,
        initial_energy_MeV=float(orbit.kinetic_MeV.iloc[0]),
        final_energy_MeV=float(orbit.kinetic_MeV.iloc[crossing]),
        final_radius_m=float(np.hypot(r[crossing, 0], r[crossing, 1])),
        first_turn_kick_steps=event_steps.tolist(),
        first_turn_kick_gains_MeV=changes[event_steps-1].tolist(),
        cavity_diagnostics_in_bounded_run=cavity_lines,
        trim_conversion='Modern named coil follows supplied cyclotron1.in; historical inline conversion not independently established.')
    (work/'summary.json').write_text(json.dumps(metrics, indent=2)+'\n')
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
