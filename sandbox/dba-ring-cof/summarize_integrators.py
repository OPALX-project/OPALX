#!/usr/bin/env python3
"""Summary of retained DBA COF grids and matched exit-search runs."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'

def main():
    rows=[]
    for method,prefix in [('BORIS','trusted-boris'),('RK4','trusted-rk4'),('DOP853','trusted-20260907')]:
        for work in sorted(OUT.glob(prefix+'-dt*')):
            p=json.loads((work/'provenance.json').read_text())
            r=json.loads((work/'closed-orbit.json').read_text())
            c=json.loads((work/'comparison.json').read_text())
            t=pd.read_csv(work/'tunes.csv').set_index('plane')
            assert r['converged']
            rows.append(dict(integrator=method,dt_s=p['dt_s'],h=p['h'],stability=r['stability'],
                qx=t.loc['horizontal','tracked_phase'],qy=t.loc['vertical','tracked_phase'],
                qx_error=t.loc['horizontal','phase_error'],qy_error=t.loc['vertical','phase_error'],
                max_tune_error=t.phase_error.abs().max(),matrix_error=c['max_matrix_error'],
                modulus_error=t.modulus_error.max(),
                position_residual=max(abs(r['residual'][i]) for i in [0,2]),
                momentum_residual=max(abs(r['residual'][i]) for i in [1,3]),
                wall_s=p['elapsed_s'],case=work.name))
    table=pd.DataFrame(rows)
    table.to_csv(OUT/'integrator-grid.csv',index=False)
    selected=table[(table.dt_s==5e-12)&(table.h==3e-5)]
    selected.to_csv(OUT/'integrator-summary.csv',index=False)
    lines=['## Integrator comparison (2026-09-07)', '',
        'Historical experiment: the Boris new-search columns use the unrestricted secant version. The final implementation restores Boris bisection; RK4/DOP853 retain secant. See README commit-scope notes and commit-boris-fallback results.', '',
        'Stable hard-edge DBA; current Release executable, one MPI rank and four OpenMP threads.',
        'All 15 COF runs verified the zero closed orbit with unchanged solver tolerances.',
        'Tables report principal fractional phases; integer and oriented branches are not determined.', '',
        '### Matched COF setting: DT=5e-12 s, h=3e-5', '',
        '| Integrator | Qx | Qy | max absolute phase error | max transverse matrix error | Classification |',
        '|---|---:|---:|---:|---:|---|']
    for row in selected.itertuples():
        lines.append(f'| {row.integrator} | {row.qx:.12f} | {row.qy:.12f} | {row.max_tune_error:.3e} | {row.matrix_error:.3e} | {row.stability} |')
    lines+=['', '### Timestep check, fixed h=3e-5', '',
            '| Integrator | DT [s] | max absolute phase error | max transverse matrix error | Classification |',
            '|---|---:|---:|---:|---|']
    for row in table[table.h==3e-5].sort_values(['integrator','dt_s'],ascending=[True,False]).itertuples():
        lines.append(f'| {row.integrator} | {row.dt_s:g} | {row.max_tune_error:.3e} | {row.matrix_error:.3e} | {row.stability} |')
    lines+=['', '### Perturbation check, fixed DT=5e-12 s', '',
            '| Integrator | h | max absolute phase error | max modulus error | Classification |',
            '|---|---:|---:|---:|---|']
    for row in table[table.dt_s==5e-12].sort_values(['integrator','h']).itertuples():
        lines.append(f'| {row.integrator} | {row.h:g} | {row.max_tune_error:.3e} | {row.modulus_error:.3e} | {row.stability} |')
    lines+=['', '### Matched old/new exit search', '',
            'DOP853, RK4 and Boris each launch their COF-verified zero orbit at the same DT=5e-12 s and h=3e-5, Richardson level 0.',
            'The direct COF comparison uses the corresponding integrator and actual beam momentum for the slope conversion.', '',
            '| Integrator | max old/new 6D difference | old/new exit trials | old/new exit seconds | old/new wall seconds |',
            '|---|---:|---:|---:|---:|']
    for method,prefix in [('BORIS','trusted-boris'),('RK4','trusted-rk4'),('DOP853','trusted-20260907')]:
        data=json.loads((OUT/(prefix+'-map-comparison.json')).read_text())
        old,new=data['rows']
        assert data['identical_design_path'] and data['identical_inputs'] and data['identical_particles']
        lines.append(f"| {method} | {data['max_old_new_difference']:.3e} | {old['exit_iterations']:.0f} / {new['exit_iterations']:.0f} | {old['ray_exit_s']:.3f} / {new['ray_exit_s']:.3f} | {old['wall_s']:.2f} / {new['wall_s']:.2f} |")
    lines+=['', 'Boris remains sensitive to exit localization: the old/new 6D difference is 1.65e-5. The new transverse map agrees more closely with the direct Boris COF Jacobian (3.13e-7 versus 1.10e-5), but the old map has a smaller transverse analytic error (7.10e-6 versus 1.68e-5). This does not establish an accuracy winner; it motivates isolating boundary/integration subdivision sensitivity before accepting the optimization for Boris.', 'Boris exhibits approximately second-order phase-error reduction over the three timesteps. Its coarse-step UNSTABLE label corresponds to only about 1.9e-8 modulus excess; it is a strict numerical diagnostic, not evidence that the analytic lattice is unstable.',
            'RK4 and DOP853 are close to the finite-difference/roundoff floor for this case; these data do not rank their asymptotic orders. Larger h adds differentiation error, and smaller h may add noise. No tolerance was relaxed.',
            'Matrix errors refer to slope coordinates, with position scaled by 1 m. Error ranges are observations against the analytic ideal-lattice reference, not rigorous uncertainty bounds. Timings are individual sequential runs. These checks do not validate the ISIS lattice or its coarse Boris step.', '',
            'Reproduce grids with run_benchmark.py --integrator METHOD --threads 4 --dt DT --h H --name FRESH_NAME; map comparisons with trusted_map_check.py --integrator METHOD --prefix FRESH_PREFIX (the named COF baseline must exist). summarize_integrators.py regenerates this report and the CSV tables from retained runs.']
    (ROOT/'INTEGRATOR_COMPARISON.md').write_text('\n'.join(lines)+'\n')
    print(selected.to_string(index=False))

if __name__=='__main__':
    main()
