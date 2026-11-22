#!/usr/bin/env python3
"""Validate the completed Lumi Track-12 output and compare with CAIN/manufactured tracks."""
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
run = HERE / 'track12'
assert (run / 'completed').exists()
assert (run / 'timing.dat').stat().st_size > 0
preparation = json.loads((run / 'results/preparation_manifest.json').read_text())
primary_hash = hashlib.sha256((run / 'input/primary_fixed.fromfile').read_bytes()).hexdigest()
assert primary_hash == preparation['primary']['sha256']
checks = []
for species in (1, 2):
    with h5py.File(run / f'track12_timed_c{species}.h5') as h5:
        previous = set()
        records = 0
        for name in sorted(h5, key=lambda item: int(item.split('#')[1])):
            group = h5[name]
            ids = group['id'][()]
            current = set(map(int, ids))
            assert len(current) == len(ids)
            assert previous <= current
            for component in ('x', 'y', 'z', 'px', 'py', 'pz'):
                assert np.isfinite(group[component][()]).all()
            previous = current
            records += len(ids)
            final_step = int(np.asarray(group.attrs['GlobalTrackStep']).item())
        assert final_step == 1500 and len(previous) == 6
        checks.append(dict(container=species, final_step=final_step,
                           final_count=len(previous), actual_h5_records=records,
                           finite=True, ids_persist=True))
pd.DataFrame(checks).to_csv(HERE / 'track12_validation.csv', index=False)
print(pd.DataFrame(checks).to_string(index=False), flush=True)
subprocess.run([sys.executable, str(ROOT / 'sandbox/experiment-2/scripts/compare_timed_track12.py'),
                '--run-dir', str(run)], check=True)
table = pd.read_csv(run / 'results/track12_pointwise_comparison.csv')
assert len(table) == 13012
report = json.loads((run / 'results/track12_comparison.json').read_text())
runtime = re.search(r'^real ([0-9.]+)$', (run / 'opalx.log').read_text(), re.M).group(1)
summary = dict(source_revision='e24dc92185198c5f65718e057400e6c6398dce56',
               slurm_job='23901513', host='lumi', mpi_ranks=1, gpu='MI250',
               primary_particles=400000, primary_sha256=primary_hash,
               mesh=[256, 256, 128], steps=1501, info_level=4,
               wall_seconds=float(runtime), validation=checks,
               comparison_rows=len(table), overall=report['overall'],
               opalx_vs_manufactured=report['opalx_vs_manufactured'])
(HERE / 'track12_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
