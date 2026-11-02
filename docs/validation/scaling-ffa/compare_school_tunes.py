#!/usr/bin/env python3
"""Validate and plot the completed production tune convergence matrix."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from school_linear_optics import optics
from school_tunes import analyze


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['coarse', 'fine', 'small', 'mpi', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    linear = optics()
    cases = {name: analyze(getattr(args, name), linear) for name in ['coarse', 'fine', 'small', 'mpi']}
    short = analyze(args.coarse, linear, window_turns=cases['mpi']['analysis_turns'])
    half = analyze(args.coarse, linear, window_turns=cases['coarse']['analysis_turns'] // 2)
    summary = dict(linear_optics=linear, cases=cases, convergence={},
                   limitation='No asymptotic second-order regime demonstrated; no zero-step extrapolation reported.')
    for plane in ['horizontal', 'vertical']:
        values = {name: value['tunes'][plane]['ring_tune'] for name, value in cases.items()}
        reference = linear[plane]['ring_tune']
        coarse_error = values['coarse'] - reference
        fine_error = values['fine'] - reference
        length_difference = abs(half['tunes'][plane]['ring_tune'] - values['coarse'])
        amplitude_difference = abs(short['tunes'][plane]['ring_tune'] - values['small'])
        mpi_difference = abs(short['tunes'][plane]['ring_tune'] - values['mpi'])
        amplitude_ratio = cases['small']['config']['amplitude_m'] / cases['coarse']['config']['amplitude_m']
        # Averaging opposite offsets removes odd terms; the leading finite-
        # amplitude tune shift is quadratic. Keep this distinct from dt^2.
        amplitude_shift = (short['tunes'][plane]['ring_tune'] - values['small']) / (1 - amplitude_ratio ** 2)
        corrected_ratio = (coarse_error - amplitude_shift) / (fine_error - amplitude_shift)
        if abs(coarse_error) >= 1e-4 or abs(fine_error) >= .4 * abs(coarse_error):
            raise AssertionError((plane, 'timestep refinement', coarse_error, fine_error))
        if length_difference >= 1e-6 or amplitude_difference >= 2e-6 or mpi_difference >= 1e-6:
            raise AssertionError((plane, 'length/amplitude/MPI', length_difference, amplitude_difference, mpi_difference))
        for case in cases.values():
            if case['tunes'][plane]['plus_minus_difference'] >= 2e-6:
                raise AssertionError((plane, 'opposite perturbations'))
        summary['convergence'][plane] = dict(coarse_linear_error=coarse_error, fine_linear_error=fine_error,
                                            error_ratio=coarse_error / fine_error,
                                            amplitude_corrected_step_error_ratio=corrected_ratio,
                                            estimated_quadratic_amplitude_shift=amplitude_shift,
                                            analysis_length_difference=length_difference,
                                            amplitude_difference=amplitude_difference, mpi_tune_difference=mpi_difference)
    # Compare actual per-ID samples across ranks for the same first complete turns.
    def common_samples(path):
        config = json.loads((path / 'config.json').read_text())
        data = pd.read_csv(path / 'probes.csv')
        ranks = config['ranks']
        data['row'] = (data.particle_id % ranks) * (8 // ranks) + data.particle_id // ranks
        data['cell'] = 16 * (data.turn - 1) + data.plane_id + np.where(data.plane_id == 0, 16, 0)
        data = data[(data.cell >= 16) & (data.cell < 16 * (cases['mpi']['analysis_turns'] + 1))]
        return data.set_index(['row', 'cell']).sort_index()
    a, b = common_samples(args.coarse), common_samples(args.mpi)
    if not a.index.equals(b.index):
        raise AssertionError('Serial/MPI sample identities differ')
    columns = ['x_m', 'y_m', 'z_m', 'px_mc', 'py_mc', 'pz_mc']
    maximum = (a[columns] - b[columns]).abs().max()
    if maximum.max() >= 1e-7:
        raise AssertionError(maximum)
    summary['serial_mpi_max_sample_difference'] = maximum.to_dict()
    for case in cases.values():
        if case['max_relative_momentum_norm_drift'] >= 1e-8:
            raise AssertionError('Momentum norm drift exceeds the diagnostic bound')
        if case['reference_duplicate_max_difference'] >= 1e-8:
            raise AssertionError('Coincident reference particles disagree')
    (args.output / 'tune-validation.json').write_text(json.dumps(summary, indent=2) + '\n')
    os.environ.setdefault('MPLCONFIGDIR', str(args.output.resolve() / '.matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(args.output.resolve() / '.cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'savefig.dpi': 220})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.1), layout='constrained')
    for ax, plane, label in zip(axes, ['horizontal', 'vertical'], ['Radial tune', 'Vertical tune']):
        target = linear[plane]['ring_tune']
        dt2 = np.array([cases[name]['config']['dt_s'] ** 2 * 1e24 for name in ['fine', 'coarse']])
        tunes = np.array([cases[name]['tunes'][plane]['ring_tune'] for name in ['fine', 'coarse']])
        ax.axhline(target, color='.4', ls='--', lw=1.2, label='Independent cell optics')
        ax.plot(dt2, tunes, '-o', color='#165b9b', label='Boris tracking, 100 µm')
        ax.scatter([1e4], [cases['small']['tunes'][plane]['ring_tune']], marker='x', color='#b95026',
                   s=50, label='10 µm perturbation', zorder=4)
        ax.set(xlabel=r'Timestep squared [ps$^2$]', ylabel='Tune per ring turn', title=label)
        ax.ticklabel_format(axis='y', style='plain', useOffset=False)
        ax.grid(alpha=.18)
    axes[0].legend(fontsize=8, loc='upper left')
    fig.suptitle('16-cell school FFA · tune convergence at 12 MeV', fontsize=12)
    for suffix in ['png', 'pdf']:
        fig.savefig(args.output / f'school-tunes.{suffix}')
    plt.close(fig)
    print(json.dumps(summary['convergence'], indent=2))


if __name__ == '__main__':
    main()
