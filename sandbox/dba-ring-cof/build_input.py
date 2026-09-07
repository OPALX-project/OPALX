#!/usr/bin/env python3
"""Generate the selected hard-edge DBA COF deck, preserving historical inputs."""
import json
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parent


def main():
    selected = json.loads((ROOT / 'analytic-results/summary.json').read_text())
    if selected['selected']['k2_m2'] is not None:
        raise ValueError('Two-family placement requires a separate geometry generator')
    builder = runpy.run_path(str(ROOT.parent / 'map-3/ring/build_ring.py'))
    text = builder['deck'](builder['geometry'](0.))
    text = text.split('DIST: DISTRIBUTION')[0]
    text = text.replace('ENABLELINEARTRANSFERMAPS=TRUE', 'ENABLELINEARTRANSFERMAPS=FALSE')
    text = text.replace('Six-cell finite-fringe DBA RING; unmatched launch',
                        'Stable non-achromatic hard-edge DBA COF benchmark')
    text = text.replace('FINT=0.1', 'FINT=0')
    text = text.replace('K1=-6.371966681365967',
                        f"K1={selected['selected']['k1_m2']:.17g}")
    text += '''BEAM1: BEAM, PARTICLE=ELECTRON, NALLOC=1,
    SOURCES="UNUSED", CHARGE=-1, PC=P0;
COF, LINE=LATTICE, BEAM=BEAM1, DT=1e-11, MAXSTEPS=20000,
    MAXPATH=35, TIMEINTEGRATOR="DOP853";
RUN, X=0, PX=0, Y=0, PY=0,
    FDSTEP={1e-5,0.0049023677597553325,1e-5,0.0049023677597553325},
    SCALES={1,490.23677597553325,1,490.23677597553325}, OUTPUT="closed-orbit";
ENDCOF;
QUIT;
'''
    output = ROOT / 'opalx'
    output.mkdir(exist_ok=True)
    (output / 'cof.in').write_text(text)
    print(output / 'cof.in')


if __name__ == '__main__':
    main()
