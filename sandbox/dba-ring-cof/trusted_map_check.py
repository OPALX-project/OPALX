#!/usr/bin/env python3
"""Matched old/new map searches on the COF-verified stable DBA orbit.

Run after run_benchmark.py's trusted DT/FD cases. Both map variants use identical
DOP853/DT/h settings and launch at the zero fixed point verified by COF.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO/'sandbox/lin-map-validation'))
from check_map import parse_map


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--integrator', choices=['BORIS','RK4','DOP853'], default='DOP853')
    parser.add_argument('--prefix', default='trusted-20260907')
    parser.add_argument('--old-binary', type=Path, default=REPO/'omp-build/src/opalx-before-exit-search')
    parser.add_argument('--new-binary', type=Path, default=REPO/'omp-build/src/opalx')
    args = parser.parse_args()
    cof_name = ('trusted-20260907-dt5e-12-h3e-5' if args.integrator == 'DOP853'
                else f'trusted-{args.integrator.lower()}-dt5e-12-h3e-5')
    cof = ROOT/'results'/cof_name
    result = json.loads((cof/'closed-orbit.json').read_text())
    assert result['converged'] and result['coordinates'] == [0,0,0,0]
    p = float(re.search(r'p0_mc=([^ ]+)', (cof/'run.log').read_text())[1])
    header = (ROOT/'opalx/cof.in').read_text().split('BEAM1: BEAM')[0]
    header = header.replace('LINEARTRANSFERMAPINTEGRATOR="DOP853"', f'LINEARTRANSFERMAPINTEGRATOR="{args.integrator}"')
    header = header.replace('ENABLELINEARTRANSFERMAPS=FALSE','ENABLELINEARTRANSFERMAPS=TRUE')
    header = header.replace('LINEARTRANSFERMAPRICHARDSON=1','LINEARTRANSFERMAPRICHARDSON=0')
    header = re.sub(r'LINEARTRANSFERMAPSTEPS=\{[^}]*\}',
                    'LINEARTRANSFERMAPSTEPS={3e-5,3e-5,3e-5,3e-5,3e-5,3e-5}',header)
    deck = header + '''DIST: DISTRIBUTION, TYPE=FROMFILE, FNAME="reference-particle.txt", NPARTDIST=1;
SOURCE: EMISSIONSOURCE, DISTRIBUTION=DIST;
SOURCES: EMISSIONSOURCELIST=(SOURCE);
FS: FIELDSOLVER, TYPE=NONE, NX=16, NY=16, NZ=16,
 PARFFTX=TRUE, PARFFTY=TRUE, PARFFTZ=TRUE, BCFFTX=OPEN, BCFFTY=OPEN, BCFFTZ=OPEN,
 BBOXINCR=1, GREENSF=INTEGRATED;
BEAM1: BEAM, PARTICLE=ELECTRON, NALLOC=1, BCHARGE=1.602176634e-19,
 SOURCES=SOURCES, CHARGE=-1;
TRACK, LINE=LATTICE, BEAM=BEAM1, MAXSTEPS=1, DT=5e-12, ZSTOP=25.766370614359175;
RUN, METHOD="PARALLEL", FIELDSOLVER=FS;
ENDTRACK;
QUIT;
'''
    analytic = np.loadtxt(ROOT/'analytic-results/ring-slope-map.txt')
    direct = np.loadtxt(cof/'tracked-slope-map.txt')
    matrices = {}
    rows = []
    for name, exe in [('old',args.old_binary),('new',args.new_binary)]:
        work = ROOT/f'results/{args.prefix}-map-{name}'
        work.mkdir(exist_ok=False)
        (work/'map.in').write_text(deck)
        particle = f'1\nx px y py z pz\n0 0 0 0 0 {p:.17g}\n'
        (work/'reference-particle.txt').write_text(particle)
        binary = exe.resolve()
        command = [str(binary),'--info','2','map.in']
        provenance = dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                          input_sha256=hashlib.sha256(deck.encode()).hexdigest(),
                          particle_sha256=hashlib.sha256(particle.encode()).hexdigest(),
                          command=command, ranks=1, omp_threads=4, p0_mc=p, integrator=args.integrator)
        (work/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
        start = time.monotonic()
        with (work/'run.log').open('w') as log:
            run = subprocess.run(command,cwd=work,stdout=log,stderr=subprocess.STDOUT,
                                 env=dict(os.environ,OMP_NUM_THREADS='4',OMP_PROC_BIND='false',OMP_DYNAMIC='false'),timeout=600)
        provenance.update(returncode=run.returncode,elapsed_s=time.monotonic()-start)
        (work/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
        assert run.returncode == 0, work
        log = (work/'run.log').read_text()
        m = np.array(parse_map(log)); matrices[name] = m
        np.savetxt(work/'slope-map.txt',m,fmt='%.16e')
        row = dict(search=name,wall_s=provenance['elapsed_s'],
                   max_6d_analytic_error=float(abs(m-analytic).max()),
                   max_transverse_analytic_error=float(abs(m[:4,:4]-analytic[:4,:4]).max()),
                   max_transverse_cof_difference=float(abs(m[:4,:4]-direct).max()))
        for label, sl in [('x',slice(0,2)),('y',slice(2,4))]:
            eig = np.linalg.eigvals(m[sl,sl])
            target = np.linalg.eigvals(analytic[sl,sl])
            row[f'q{label}'] = float(max(abs(np.angle(eig)))/(2*np.pi))
            row[f'q{label}_error'] = row[f'q{label}']-float(max(abs(np.angle(target)))/(2*np.pi))
            row[f'{label}_modulus_error'] = float(max(abs(abs(eig)-1)))
        for key in ['exit_iterations','ray_exit_s','field_samples','support_lookups']:
            matches = re.findall(r'\* OT maps: .*?\b'+key+r'=([\d.e+-]+)',log)
            row[key] = float(matches[-1])
        rows.append(row)
        print(json.dumps(row,indent=2),flush=True)
    summary = dict(rows=rows,max_old_new_difference=float(abs(matrices['old']-matrices['new']).max()))
    old = ROOT/f'results/{args.prefix}-map-old'
    new = ROOT/f'results/{args.prefix}-map-new'
    a = json.loads((old/'provenance.json').read_text())
    b = json.loads((new/'provenance.json').read_text())
    assert a['input_sha256'] == b['input_sha256']
    assert a['particle_sha256'] == b['particle_sha256']
    paths = [pd.read_csv(next((w/'data').glob('*DesignPath.dat')),sep=r'\s+',
                         comment='#',header=None,usecols=range(15)) for w in [old,new]]
    assert paths[0].equals(paths[1])
    pattern = r'^OPAL-X> \* RING (?:reference return length|return displacement).*'
    assert re.findall(pattern,(old/'run.log').read_text(),re.M) == re.findall(pattern,(new/'run.log').read_text(),re.M)
    summary.update(identical_inputs=True,identical_particles=True,identical_design_path=True,
                   identical_printed_reference_return=True,design_path_rows=len(paths[0]))
    (ROOT/f'results/{args.prefix}-map-comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame(rows).to_csv(ROOT/f'results/{args.prefix}-map-comparison.csv',index=False)

if __name__ == '__main__':
    main()
