#!/usr/bin/env python3
"""Collect the retained single-rank trusted DBA study without rerunning tracking."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
from analytic_scan import cell

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT/'results'

def main():
    # Verify that the retained reference is reproduced by the analytic formulas.
    analytic = np.linalg.matrix_power(cell(2.84), 6)
    np.testing.assert_allclose(np.loadtxt(ROOT/'analytic-results/ring-slope-map.txt'),
                               analytic, atol=1e-12, rtol=1e-13)
    rows = []
    for work in sorted(RESULTS.glob('trusted-20260907-dt*')):
        result = json.loads((work/'comparison.json').read_text())
        provenance = json.loads((work/'provenance.json').read_text())
        tunes = pd.read_csv(work/'tunes.csv').set_index('plane')
        rows.append(dict(case=work.name,dt_s=provenance['dt_s'],h=provenance['h'],
                         stability=result['stability'],matrix_error=result['max_matrix_error'],
                         max_position_residual_m=max(abs(result['residual'][i]) for i in [0,2]),
                         max_momentum_residual=max(abs(result['residual'][i]) for i in [1,3]),
                         qx=tunes.loc['horizontal','tracked_phase'],qy=tunes.loc['vertical','tracked_phase'],
                         qx_error=tunes.loc['horizontal','phase_error'],qy_error=tunes.loc['vertical','phase_error'],
                         max_modulus_error=tunes.modulus_error.max()))
    table = pd.DataFrame(rows).sort_values(['h','dt_s'])
    table.to_csv(RESULTS/'trusted-20260907-cof-summary.csv',index=False)
    print(table.to_string(index=False))
    selected = table[(table.dt_s==5e-12)&(table.h==3e-5)].iloc[0].to_dict()
    summary = dict(selected=selected,
                   observed_max_tune_error=float(table[['qx_error','qy_error']].abs().max().max()),
                   observed_max_matrix_error=float(table.matrix_error.max()),
                   note='Observed envelope against the analytic reference over the tested settings; not a rigorous error bound or an asymptotic convergence-order measurement.')
    (RESULTS/'trusted-20260907-cof-summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__ == '__main__':
    main()
