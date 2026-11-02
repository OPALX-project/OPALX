#!/usr/bin/env python3
"""Check saved COF runs against an independent polar ODE and track returns."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from school_reference import solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reference, _ = solve()
    expected = np.array([reference['x_m'], reference['radial_momentum_mc'], 0., 0.])
    reports = []
    for case in args.cases:
        result = json.loads((case / 'orbit.json').read_text())
        log = (case / 'run.log').read_text()
        assert result['converged'] and abs(result['energy_MeV'] - 12.) < 1e-9
        u, residual = np.array(result['coordinates']), np.array(result['residual'])
        assert np.max(np.abs(residual[[0, 2]])) < 1e-10
        assert np.max(np.abs(residual[[1, 3]])) < 1e-11
        error = u - expected
        # 1 nm and 1e-10 mc include the configured Newton tolerance and
        # integration differences; the independent ODE closes below 1 pm.
        assert np.max(np.abs(error[[0, 2]])) < 1e-9, error
        assert np.max(np.abs(error[[1, 3]])) < 1e-10, error
        assert abs(result['relative_energy_drift']) < 5e-9
        evaluations = re.search(r'evaluations=(\d+)', log)
        returned = re.search(r'Return time_s=(\S+) path_m=(\S+)', log)
        integrator = re.search(r'Integrator=(\S+)', log)[1]
        ranks = int(re.search(r'MPI ranks: (\d+)', log)[1])
        fdsteps = [float(v) for v in re.search(r'FDSTEP=([^\n]+)', log)[1].split()]
        record = dict(case=case.name, dt_s=result['dt_s'], integrator=integrator, ranks=ranks,
                      finite_difference_steps=fdsteps, coordinates=u.tolist(),
                      residual=residual.tolist(), reference_error=error.tolist(),
                      relative_energy_drift=result['relative_energy_drift'],
                      evaluations=int(evaluations[1]), return_time_s=float(returned[1]),
                      return_path_m=float(returned[2]))
        if (case / 'probes.csv').exists():
            assert 'completed directed turn 3,' in log
            data = pd.read_csv(case / 'probes.csv')
            data = data[(data.status == 'crossing') & (data.turn <= 2)]
            assert not data.duplicated(['particle_id', 'plane_id', 'turn']).any()
            counts = data.groupby('particle_id').size()
            assert len(counts) == 8 and (counts == 2).all()
            position_error = float((data.x_m - u[0]).abs().max())
            momentum_error = float((data.px_mc - u[1]).abs().max())
            assert position_error < 5e-6 and momentum_error < 5e-7
            record['tracking'] = dict(complete_compared_turns=2, particles=8,
                                      max_return_x_error_m=position_error,
                                      max_return_px_error_mc=momentum_error)
        reports.append(record)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary = dict(reference=reference, runs=reports)
    args.output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
