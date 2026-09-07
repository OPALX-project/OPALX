#!/usr/bin/env python3
"""Single-rank DBA COF comparison; fresh directories preserve prior evidence."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True)
    parser.add_argument('--dt', type=float, default=1e-11)
    parser.add_argument('--h', type=float, default=1e-5, help='position m / slope FD amplitude')
    parser.add_argument('--integrator', choices=['BORIS','RK4','DOP853'], default='DOP853')
    parser.add_argument('--displaced', action='store_true')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--binary', type=Path, default=REPO / 'omp-build/src/opalx')
    args = parser.parse_args()
    work = ROOT / 'results' / args.name
    work.mkdir(parents=True, exist_ok=False)
    text = (ROOT / 'opalx/cof.in').read_text()
    text = text.replace('DT=1e-11', f'DT={args.dt:.17g}')
    text = text.replace('TIMEINTEGRATOR="DOP853"', f'TIMEINTEGRATOR="{args.integrator}"')
    p = 490.23677663749817
    text = re.sub(r'FDSTEP=\{[^}]*\}',
                  f'FDSTEP={{{args.h},{p*args.h},{args.h},{p*args.h}}}', text)
    if args.displaced:
        text = text.replace('RUN, X=0, PX=0, Y=0, PY=0,',
                            'RUN, X=0.0001, PX=0.001, Y=0.0001, PY=0.001,')
    (work / 'cof.in').write_text(text)
    binary = args.binary.resolve()
    provenance = dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                      input_sha256=hashlib.sha256(text.encode()).hexdigest(), ranks=1, omp_threads=args.threads, dt_s=args.dt, h=args.h, integrator=args.integrator)
    start = time.monotonic()
    with (work / 'run.log').open('w') as log:
        run = subprocess.run([str(binary), '--info', '1', 'cof.in'], cwd=work,
                             env=dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_PROC_BIND='false', OMP_DYNAMIC='false'),
                             stdout=log, stderr=subprocess.STDOUT, timeout=600)
    provenance.update(returncode=run.returncode, elapsed_s=time.monotonic()-start)
    (work / 'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    if run.returncode:
        raise RuntimeError((work / 'run.log').read_text()[-5000:])
    result = json.loads((work / 'closed-orbit.json').read_text())
    p = float(re.search(r'p0_mc=([^ ]+)', (work / 'run.log').read_text())[1])
    transform = np.diag([1., p, 1., p])
    measured = np.linalg.inv(transform) @ np.array(result['matrix']) @ transform
    analytic = np.loadtxt(ROOT / 'analytic-results/ring-slope-map.txt')[:4, :4]
    np.savetxt(work / 'tracked-slope-map.txt', measured, fmt='%.16e')
    np.savetxt(work / 'matrix-difference.txt', measured-analytic, fmt='%.16e')
    rows = []
    for label, block in [('horizontal', slice(0,2)), ('vertical', slice(2,4))]:
        target = np.linalg.eigvals(analytic[block, block])
        actual = np.linalg.eigvals(measured[block, block])
        q = max(abs(np.angle(actual)))/(2*np.pi)
        expected = max(abs(np.angle(target)))/(2*np.pi)
        rows.append(dict(plane=label, analytic_phase=expected, tracked_phase=q,
                         phase_error=q-expected,
                         modulus_error=max(abs(abs(actual)-1))))
    table = pd.DataFrame(rows)
    table.to_csv(work / 'tunes.csv', index=False)
    summary = dict(max_matrix_error=float(abs(measured-analytic).max()),
                   coordinates=result['coordinates'], residual=result['residual'],
                   stability=result['stability'], relative_energy_drift=result['relative_energy_drift'])
    (work / 'comparison.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(table.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
