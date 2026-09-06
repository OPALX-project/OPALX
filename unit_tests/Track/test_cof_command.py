#!/usr/bin/env python3
"""End-to-end COF parser/runtime regressions using native magnetic ring inputs.

Only the Python standard library is needed. Temporary artifacts are isolated.
MPI comparisons are dedicated regressions, not physics convergence studies.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--opalx', type=Path, required=True)
parser.add_argument('--mpiexec')
args = parser.parse_args()
exe = str(args.opalx.resolve())
mpi_flags = []
if args.mpiexec and 'Open MPI' in subprocess.check_output([args.mpiexec, '--version'], text=True):
    mpi_flags = ['--map-by', 'slot', '--bind-to', 'none', '--oversubscribe']
base = '''OPTION, VERSION=10900;
REAL P0=0.5;
B: BEAM, PARTICLE=PROTON, PC=P0, NALLOC=1, BCHARGE=1.602176634e-19, SOURCES="UNUSED";
D: DRIFT, L=1, ELEMEDGE=0;
BAD: LINE=(D);
B0: SBEND, L=PI/2, ANGLE=PI/2, K1=-0.3, ELEMEDGE=0;
B1: SBEND, L=PI/2, ANGLE=PI/2, K1=-0.3, ELEMEDGE=PI/2;
B2: SBEND, L=PI/2, ANGLE=PI/2, K1=-0.3, ELEMEDGE=PI;
B3: SBEND, L=PI/2, ANGLE=PI/2, K1=-0.3, ELEMEDGE=3*PI/2;
R: RING=(B0,B1,B2,B3);
'''
block = '''COF, LINE=R, BEAM=B, DT=1e-11, MAXSTEPS=10000, MAXPATH=10;
RUN, X=0.0001, Y=0.0001, OUTPUT="result";
ENDCOF;
'''
with tempfile.TemporaryDirectory(prefix='opalx-cof-command-') as tmp:
    root = Path(tmp)

    def run(name, text, expected=None, ranks=1):
        work = root / name
        work.mkdir()
        (work / 'test.in').write_text(text)
        command = [exe, '--info', '1', 'test.in']
        if ranks > 1:
            command = [args.mpiexec] + mpi_flags + ['-n', str(ranks)] + command
        result = subprocess.run(command, cwd=work, text=True,
                                env=dict(os.environ, OMP_NUM_THREADS='1', OMP_PROC_BIND='false',
                                         OMPI_MCA_rmaps_base_oversubscribe='1',
                                         OMPI_MCA_rmaps_default_mapping_policy='slot'),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
        diagnostics = result.stdout + result.stderr
        if expected:
            assert expected in diagnostics, (name, diagnostics)
        else:
            assert result.returncode == 0 and (work / 'result.json').exists(), (name, diagnostics)
            # Check stdout separately from stderr and files, including MPI runs.
            # Exactly one report per block: worker ranks must not duplicate it.
            blocks = text.count('COF,')
            for heading in ('Coordinates (x,px,y,py):', 'Matrix:',
                            'One-turn eigenanalysis:', 'lambda[0]', 'lambda[1]',
                            'lambda[2]', 'lambda[3]', 'Mode 1:'):
                assert result.stdout.count(heading) == blocks, (name, heading, diagnostics)
            value = json.loads((work / 'result.json').read_text())
            assert value['converged'] and max(map(abs, value['residual'])) <= 1e-10
            assert value['energy_MeV'] > 0
            assert len(value['eigenvalues']) == 4
            assert len((work / 'result-orbit.csv').read_text().splitlines()) > 100
        print(name, 'PASS', flush=True)
        return work

    one = run('native-ring', base + block)
    run('electron-ring', base.replace('PARTICLE=PROTON', 'PARTICLE=ELECTRON') + block)
    run('aperture-loss', base.replace('K1=-0.3,', 'K1=-0.3, HAPERT=0.001, HGAP=0.01,')
        + block.replace('X=0.0001', 'X=0.1'), 'Aperture/material loss')
    # Opposite drifts share longitudinal slabs but are distinct ring arms.
    # Never let the remote slab reject the local orbit; still reject a launch
    # outside the aperture of its nearest finite design segment.
    racetrack = base.split('D: DRIFT')[0]
    members = []
    for i in range(4):
        start = i*(math.pi/2+1)
        racetrack += (f'RB{i}: SBEND, L=PI/2, ANGLE=PI/2, K1=-0.3, ELEMEDGE={start};\n'
                      f'RD{i}: DRIFT, L=1, APERTURE="RECTANGLE(0.002,0.002)", '
                      f'ELEMEDGE={start+math.pi/2};\n')
        members += [f'RB{i}', f'RD{i}']
    racetrack += 'R: RING=('+','.join(members)+');\n'
    drift_block = block.replace('MAXPATH=10',
        'MAXPATH=12, SECTION={-1.5,0,1,-PI/2,0,0}')
    run('opposite-drift-arms', racetrack + drift_block)
    run('local-drift-aperture-loss', racetrack + drift_block.replace('X=0.0001', 'X=0.01'),
        'Aperture/material loss in RD0')
    run('line-rejected', base + block.replace('LINE=R', 'LINE=BAD'), 'LINE must name a RING')
    run('geometry-rejected', base.replace('RING=(B0,B1,B2,B3)', 'RING=(B0,B1,B2)') + block,
        'Nominal ring geometry does not close')
    run('missing-energy', base.replace('PC=P0, ', '') + block, 'BEAM requires explicit')
    run('invalid-dt', base + block.replace('DT=1e-11', 'DT=-1'), 'DT must be finite and positive')
    run('fractional-steps', base + block.replace('MAXSTEPS=10000', 'MAXSTEPS=1.5'), 'positive integer')
    run('bad-fd', base + block.replace('RUN,', 'RUN, FDSTEP={1e-5,1e-5},'), 'exactly four')
    run('invalid-pz', base + block.replace('X=0.0001', 'PX=10'), 'Launch has no real forward')
    run('no-return', base + block.replace('MAXSTEPS=10000', 'MAXSTEPS=1'), 'No directed return')
    run('missing-end', base + block.replace('ENDCOF;', ''), 'COF requires one RUN followed by ENDCOF')
    run('empty-block', base + block.replace('RUN, X=0.0001, Y=0.0001, OUTPUT="result";', ''),
        'COF requires one RUN followed by ENDCOF')
    run('rf-rejected', base.replace('R: RING=', 'RF: RFCAVITY, L=0.1, VOLT=0, FREQ=50, ELEMEDGE=0;\nR: RING=').replace(
        'RING=(B0,B1,B2,B3)', 'RING=(B0,B1,B2,B3,RF)') + block, 'Unsupported element')
    two_blocks = block + block.replace('OUTPUT="result"', 'OUTPUT="second"')
    scoped = run('scoped-blocks', base + two_blocks)
    assert (scoped / 'result.json').read_bytes() == (scoped / 'second.json').read_bytes()
    run('output-protected', base + block + block, 'COF output exists')
    run('explicit-section', base.replace('RING=(B0,B1,B2,B3);',
        'RING=(B0,B1,B2,B3), X=2, Y=3, Z=4, THETA=0.2;')
        + block.replace('MAXPATH=10', 'MAXPATH=10, SECTION={2,3,4,0.2,0,0}'))
    if args.mpiexec:
        two = run('mpi-two', base + block, ranks=2)
        assert (one / 'result.json').read_bytes() == (two / 'result.json').read_bytes()
        run('mpi-failure', base + block.replace('MAXSTEPS=10000', 'MAXSTEPS=1'),
            'No directed return', ranks=2)
