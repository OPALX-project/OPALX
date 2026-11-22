#!/usr/bin/env python3
"""Record actual Lumi provenance and reanalyze the fetched smoke outputs locally."""
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
fingerprint = (HERE / 'execution-fingerprint.txt').read_text()
revision = fingerprint.splitlines()[0].strip()
exe_hash = re.search(r'^([0-9a-f]{64})\s+.*?/opalx$', fingerprint, re.M).group(1)
job = re.search(r'JobId=(\d+)', fingerprint).group(1)
rows = []
for experiment in range(1, 5):
    run = ROOT / f'sandbox/experiment-{experiment}/results/lumi-e24dc9218-smoke'
    manifest = run / 'run_manifest.json'
    preparation = run / 'preparation_run_manifest.json'
    if not preparation.exists():
        preparation.write_bytes(manifest.read_bytes())
    data = json.loads(preparation.read_text())
    data['preparation_executable_sha256'] = data.pop('executable_sha256', None)
    data['executable_sha256'] = exe_hash
    data['git_revision'] = revision
    data['git_status'] = ''
    data['remote_execution'] = {
        'host': 'lumi', 'slurm_job': job, 'backend': 'HIP/GFX90A',
        'mpi_ranks': 1, 'MPICH_GPU_IPC_ENABLED': '0',
        'run_directory': f'/flash/project_465002453/bb-validation-e24dc9218-20261120/experiment-{experiment}',
        'fingerprint': str(HERE / 'execution-fingerprint.txt'),
    }
    manifest.write_text(json.dumps(data, indent=2) + '\n')
    subprocess.run([sys.executable, str(ROOT / 'sandbox/scripts/run_experiment.py'),
                    '--experiment', str(experiment), '--analyze-only', str(run)], check=True)
    result = json.loads(manifest.read_text())
    rows.append(dict(experiment=experiment, status=result['status'], slurm_job=job))
    if experiment == 1:
        table = pd.read_csv(run / 'analysis/opalx_vs_rigid_gaussian_summary.csv')
        assert (table.uncovered_probe_fraction == 0).all()
        print(table[['case', 'relative_l2_E', 'relative_l2_B', 'uncovered_probe_fraction']])
    if experiment == 3:
        for species in (1, 2):
            with h5py.File(run / f'gamma_gamma_pairs_timed_c{species}.h5') as h5:
                first = h5['Step#0']
                momentum = np.array([first[k][0] for k in ('px', 'py', 'pz')])
                print('First witness momentum [mc]', species, momentum,
                      'magnitude', np.linalg.norm(momentum))
                previous = set()
                for step in sorted(h5, key=lambda name: int(name.split('#')[1])):
                    values = h5[step]['id'][()]
                    current = set(map(int, values))
                    assert len(current) == len(values), 'Duplicate IDs'
                    assert previous <= current, 'Previously emitted ID disappeared'
                    previous = current
                assert len(previous) == 1297
pd.DataFrame(rows).to_csv(HERE / 'smoke_status.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
