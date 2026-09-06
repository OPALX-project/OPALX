#!/usr/bin/env python3
"""Validate and plot the completed one-rank 590 MeV code comparisons.

Acceptance is inter-code agreement, not equality to exactly 590 MeV or an
absolute-accuracy claim. The 1.5 keV / 100 um envelopes cover measured legacy
gap-time-estimate differences; refinement must reduce the no-TC discrepancy.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    root = Path(__file__).resolve().parent/'acceleration590'
    cases = [('No TC, DT', 'opalx-endpoint-2880-no-tc'),
             ('No TC, DT/2', 'opalx-endpoint-5760-no-tc'),
             ('TC on, DT', 'opalx-endpoint-2880-tc')]
    rows = []
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    colors = ['#0072B2', '#D55E00', '#009E73']
    for (name, directory), color in zip(cases, colors):
        work = root/directory
        result = json.loads((work/'summary.json').read_text())
        assert result['ranks'] == 1
        assert result['final_energy_MeV'] >= 590
        assert abs(result['final_energy_difference_eV']) < 1500, result
        assert result['max_sampled_position_error_m'] < 1e-4, result
        assert result['final_momentum_error'] < 3e-6, result
        assert result['max_particle_offset_m'] < 1e-7, result
        assert result['gap_event_count'] == 949, result
        rows.append(dict(case=name, old_MeV=result['old_final_energy_MeV'],
            opalx_MeV=result['final_energy_MeV'], difference_eV=result['final_energy_difference_eV'],
            difference_ppm=result['final_energy_difference_eV']/result['old_final_energy_MeV'],
            final_position_error_um=result['final_position_error_m']*1e6,
            max_position_error_um=result['max_sampled_position_error_m']*1e6,
            particle_reference_nm=result['max_particle_offset_m']*1e9))
        energy = pd.read_csv(work/'energy.csv')
        error = pd.read_csv(work/'errors.csv')
        axes[0].plot(energy.time_s*1e6, energy.energy_MeV, color=color, lw=1.2, label=name)
        axes[1].plot(error.time_s*1e6, error.position_error_m*1e6, color=color, lw=1.2, label=name)
    table = pd.DataFrame(rows)
    assert abs(table.difference_eV.iloc[1]) < abs(table.difference_eV.iloc[0])
    assert table.max_position_error_um.iloc[1] < table.max_position_error_um.iloc[0]
    assert abs(table.old_MeV.iloc[2]-table.old_MeV.iloc[0]) > .01, 'TC-on must have a measurable effect'
    table.to_csv(root/'comparison.csv', index=False)
    (root/'comparison.json').write_text(json.dumps(rows, indent=2)+'\n')
    axes[0].set(xlabel='Time [µs]', ylabel='Kinetic energy [MeV]', title='72 to 590 MeV acceleration')
    axes[1].set(xlabel='Time [µs]', ylabel='Position difference [µm]', title='OPALX minus old OPAL LF2')
    axes[0].axhline(590, color='0.5', lw=.7, ls='--')
    for ax in axes: ax.legend(frameon=False, fontsize=8)
    fig.savefig(root/'comparison.png', dpi=200)
    fig.savefig(root/'comparison.pdf')
    print(table.to_string(index=False))


if __name__ == '__main__':
    main()
