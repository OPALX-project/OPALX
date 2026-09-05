#!/usr/bin/env python3
"""Generate and run one-rank OPALX acceleration; compare to the old LF2 first turn."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps-per-turn', type=int, default=2880)
    parser.add_argument('--no-trim-coils', action='store_true')
    parser.add_argument('--target-mev', type=float)
    parser.add_argument('--max-turns', type=int, default=220)
    parser.add_argument('--fixed-steps', type=int, help='Compare a matched endpoint; do not enable EKINSTOP')
    parser.add_argument('--reference', type=Path,
        default=root/'acceleration-reference/lf2-2880-validated/first-turn.csv')
    args = parser.parse_args()
    work = args.output.resolve()
    work.mkdir(parents=True, exist_ok=False)
    old = pd.read_csv(args.reference)
    launch = old.iloc[0]
    x0, z0 = launch.x, launch.y
    lines = ['OPTION, VERSION=10900, PSDUMPFREQ=1, STATDUMPFREQ=1;',
        'OPTION, ENABLELINEARTRANSFERMAPS=FALSE;',
        'TC1: TRIMCOIL, TYPE="PSI-BFIELD-MIRRORED", RMIN=4.350, RMAX=4.470, BMAX=1.4e-3, SLPTC=600;']
    if args.no_trim_coils:
        lines[-1] = lines[-1].replace('BMAX=1.4e-3', 'BMAX=0')
    if args.target_mev:
        lines[0] = 'OPTION, VERSION=10900, PSDUMPFREQ=1000, STATDUMPFREQ=1000, STEPINFOFQ=10000;'
    for i in range(8):
        a = i*math.pi/4
        lines.append(f'SM{i}: CYCLOTRONSECTOR, FMAPFN="{root / "opal/bfield.dat"}", '
            f'SYM=8, TRIMCOIL={{"TC1"}}, X={3.3*math.cos(a)-x0:.17g}, '
            f'Z={3.3*math.sin(a)-z0:.17g}, THETA={-a:.17g};')
    for name, angle, voltage, freq, phase, lo, hi, offset, width, filename in [
        ('RF0',35,.847,50.65,204,1.00338,6.20338,.416,.3,'rffield1.dat'),
        ('RF1',125,.847,50.65,384,1.00338,6.20338,.416,.3,'rffield1.dat'),
        ('RF3',215,.847,50.65,204,1.00338,6.20338,.416,.3,'rffield1.dat'),
        ('RF4',305,.847,50.65,384,1.00338,6.20338,.416,.3,'rffield1.dat'),
        ('RF2',260,.847*4*.112*1.465,151.95,149.4,1.84545,4.47724,.452,.25,'rffield2.dat')]:
        a = math.radians(angle)
        # Old plane x*sin(a)-y*cos(a)=PDIS: origin lies on -PDIS*tangent.
        lines.append(f'{name}: RFCAVITY, TYPE="SINGLEGAP", FMAPFN="{root / "opal" / filename}", '
            f'VOLT={voltage:.17g}, FREQ={freq}, PHI0={phase}, RMIN={lo}, RMAX={hi}, '
            f'GAPWIDTH={width}, X={offset*math.sin(a)-x0:.17g}, '
            f'Z={-offset*math.cos(a)-z0:.17g}, THETA={-a:.17g};')
    names = ','.join([f'SM{i}' for i in range(8)]+['RF0','RF1','RF2','RF3','RF4'])
    lines += [f'MYCYCL: RING=({names}), X={x0:.17g}, Z={z0:.17g};',
        'DIST0: DISTRIBUTION, TYPE=FROMFILE, FNAME="proton.dat", NPARTDIST=1;',
        'ES0: EMISSIONSOURCE, DISTRIBUTION=DIST0;', 'SOURCES0: EMISSIONSOURCELIST=(ES0);',
        'FS0: FIELDSOLVER, TYPE=NONE, NX=16, NY=16, NZ=16, PARFFTX=TRUE, PARFFTY=TRUE, PARFFTZ=TRUE, BCFFTX=OPEN, BCFFTY=OPEN, BCFFTZ=OPEN;',
        'BEAM0: BEAM, PARTICLE=PROTON, NALLOC=1, BCHARGE=1.602176634e-19, SOURCES=SOURCES0, CHARGE=1;',
        f'TRACK, LINE=MYCYCL, BEAM=BEAM0, MAXSTEPS={args.fixed_steps or args.steps_per_turn*(args.max_turns if args.target_mev else 2)}, DT={6/(50.65e6*args.steps_per_turn):.17g}, ZSTOP=1e6'
        + (f', EKINSTOP={args.target_mev/1000:.17g}' if args.target_mev and not args.fixed_steps else '') + ';',
        'RUN, METHOD="PARALLEL", FIELDSOLVER=FS0' + ('' if args.target_mev or args.fixed_steps else ', TURNS=1') + ';', 'ENDTRACK;', 'QUIT;']
    (work/'acceleration.in').write_text('\n'.join(lines)+'\n')
    (work/'proton.dat').write_text(f'1\nx px y py z pz\n0 {launch.px:.17g} 0 0 0 {launch.py:.17g}\n')
    with (work/'run.log').open('w') as log:
        subprocess.run([str(root.parents[1]/'omp-build/src/opalx'), '--info', '2', 'acceleration.in'],
            cwd=work, env=dict(os.environ, OMP_NUM_THREADS='1'), stdout=log,
            stderr=subprocess.STDOUT, check=True)
    rows = []
    with h5py.File(work/'acceleration.h5') as f:
        for g in f.values():
            if 'RefPartR' not in g.attrs: continue
            rows.append(dict(t=float(g.attrs['TIME'][0]), energy=float(g.attrs['ENERGY'][0]),
                r=np.array(g.attrs['RefPartR']), p=np.array(g.attrs['RefPartP']),
                offset=float(np.linalg.norm([g[k][0] for k in ['x','y','z']]))))
    rows.sort(key=lambda r:r['t'])
    events = [s for s in (work/'run.log').read_text().splitlines() if 'Cyclotron gap ' in s]
    result = dict(steps=len(rows), final_energy_MeV=rows[-1]['energy'],
                  final_time_s=rows[-1]['t'], max_particle_offset_m=max(r['offset'] for r in rows),
                  gap_events=events)
    if args.target_mev:
        result['output_records'] = result.pop('steps')
        result['position_m'] = rows[-1]['r'].tolist()
        result['momentum'] = rows[-1]['p'].tolist()
        result['target_MeV'] = args.target_mev
        result['ranks'] = 1
        result['trim_coils'] = not args.no_trim_coils
        if not args.fixed_steps:
            assert 'EKINSTOP reached after complete gap kick' in (work/'run.log').read_text()
        assert result['final_energy_MeV'] >= args.target_mev
        assert result['max_particle_offset_m'] < 1e-7, result
        result['gap_event_count'] = len(result.pop('gap_events'))
        if args.fixed_steps:
            indices = np.rint(np.array([r['t'] for r in rows])*50.65e6*args.steps_per_turn/6).astype(int)
            matched = old.iloc[indices]
            errors_r = np.linalg.norm(np.array([r['r'] for r in rows])-matched[['x','z','y']].to_numpy()*[1,-1,1], axis=1)
            errors_p = np.linalg.norm(np.array([r['p'] for r in rows])-matched[['px','pz','py']].to_numpy()*[1,-1,1], axis=1)
            result['max_sampled_position_error_m'] = float(errors_r.max())
            result['final_position_error_m'] = float(errors_r[-1])
            result['final_momentum_error'] = float(errors_p[-1])
            result['old_final_energy_MeV'] = float(matched.kinetic_MeV.iloc[-1])
            result['final_energy_difference_eV'] = (result['final_energy_MeV']-result['old_final_energy_MeV'])*1e6
            pd.DataFrame({'time_s':[r['t'] for r in rows], 'position_error_m':errors_r, 'momentum_error':errors_p}).to_csv(work/'errors.csv',index=False)
        (work/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
        pd.DataFrame({'time_s':[r['t'] for r in rows], 'energy_MeV':[r['energy'] for r in rows]}).to_csv(work/'energy.csv',index=False)
        print(json.dumps(result,indent=2))
        return
    n = min(len(rows), len(old)-1)
    reference_r = old[['x','z','y']].to_numpy()[1:n+1]*[1,-1,1]
    reference_p = old[['px','pz','py']].to_numpy()[1:n+1]*[1,-1,1]
    times = np.array([r['t'] for r in rows[:n]])
    if not np.allclose(times, old.time_s.to_numpy()[1:n+1], rtol=1e-10, atol=1e-18):
        raise RuntimeError('Use an old reference with the same timestep')
    errors = pd.DataFrame(dict(time_s=times,
        position_error_m=np.linalg.norm(np.array([r['r'] for r in rows[:n]])-reference_r,axis=1),
        momentum_error=np.linalg.norm(np.array([r['p'] for r in rows[:n]])-reference_p,axis=1)))
    errors.to_csv(work/'errors.csv', index=False)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    new_r = np.array([r['r'] for r in rows])
    axes[0].plot(old.x, old.y, color='0.5', lw=2, label='OPAL 2022.1 LF2')
    axes[0].plot(new_r[:,0], new_r[:,2], color='#0072B2', lw=1, ls='--', label='OPALX')
    axes[0].set(xlabel='X [m]', ylabel='Z [m]', title='First accelerating turn', aspect='equal')
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(times*1e9, errors.position_error_m*1e6, color='#D55E00')
    axes[1].set(xlabel='Time [ns]', ylabel='Position difference [µm]',
                title='Matched-time orbit comparison')
    fig.savefig(work/'comparison.png', dpi=200)
    fig.savefig(work/'comparison.pdf')
    plt.close(fig)
    result['max_position_error_m'] = float(errors.position_error_m.max())
    result['max_momentum_error'] = float(errors.momentum_error.max())
    result['old_final_energy_MeV'] = float(old.kinetic_MeV.iloc[-1])
    assert 'completed directed turn 1' in (work/'run.log').read_text()
    assert len(events) == 5, events
    assert result['max_particle_offset_m'] < 1e-8, result
    # First-turn compatibility envelope: the localized crossing intentionally
    # differs from old OPAL's distance/speed estimate (observed 0.76 um at DT).
    # These are regression guards, not an asymptotic convergence claim.
    assert result['max_position_error_m'] < 2e-6, result
    assert result['max_momentum_error'] < 5e-7, result
    assert abs(result['final_energy_MeV']-result['old_final_energy_MeV']) < 2e-5, result
    (work/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
