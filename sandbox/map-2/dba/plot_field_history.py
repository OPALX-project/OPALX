#!/usr/bin/env python3
"""Record and plot the DBA reference magnetic field from an ASCII SDDS .stat file.

The convergence grid uses MAXSTEPS=1 and cannot provide a time history in .stat.
This separate diagnostic copy extends tracking and records statistics every step;
it is not a convergence-grid point. Fields in .stat belong to the production
reference trajectory (Boris), not the private RK map rays.
"""
from __future__ import annotations

import argparse
import io
import json
import os
from pathlib import Path
import re
import subprocess

from convergence_dt import np, pd, plt
from convergence_grid import digest
from convergence_perturbations import configure_map


def read_ascii_sdds(path: Path) -> pd.DataFrame:
    """Read the no-row-count ASCII SDDS format emitted by OPALX writers."""
    text = path.read_text()
    blocks = re.findall(r'&(column|parameter|data)\b(.*?)&end', text, flags=re.S | re.I)
    columns, units, parameters = [], {}, 0
    for kind, body in blocks:
        if kind.lower() == 'parameter':
            parameters += 1
        elif kind.lower() == 'column':
            name = re.search(r'\bname\s*=\s*"?([^",\s]+)', body).group(1)
            unit = re.search(r'\bunits\s*=\s*"?([^",\s]+)', body)
            columns.append(name)
            units[name] = unit.group(1) if unit else None
        else:
            if not re.search(r'\bmode\s*=\s*ascii', body, re.I):
                raise ValueError('expected ASCII SDDS')
            if not re.search(r'\bno_row_counts\s*=\s*1', body):
                raise ValueError('expected no_row_counts=1')
    marker = re.search(r'&data\b.*?&end', text, flags=re.S | re.I)
    if marker is None:
        raise ValueError('missing SDDS data declaration')
    rows = [line for line in text[marker.end():].splitlines() if line.strip() and not line.lstrip().startswith('!')]
    frame = pd.read_csv(io.StringIO('\n'.join(rows[parameters:])), sep=r'\s+', names=columns,
                        keep_default_na=False)
    frame.attrs['units'] = units
    return frame


def read_stat(path: Path) -> pd.DataFrame:
    frame = read_ascii_sdds(path)
    units = frame.attrs['units']
    expected = {'t': 'ns', 's': 'm', 'Bx_ref': 'T', 'By_ref': 'T', 'Bz_ref': 'T'}
    if any(units.get(key) != unit for key, unit in expected.items()):
        raise ValueError(f'unexpected units: {units}')
    if (len(frame) < 2 or not frame.t.is_monotonic_increasing or not frame.s.is_monotonic_increasing
            or not np.isfinite(frame[list(expected)]).all().all()):
        raise ValueError('need a finite, increasing field time history with at least two samples')
    return frame


def read_element_spans(path: Path) -> pd.DataFrame:
    """Decode IndexMap step rows into separate, possibly overlapping spans.

    The last row at a repeated s gives the names active on the following open
    interval. Zero-width rows are transitions, not finite elements. Adjacent
    intervals for one name are merged; separated visits are never bridged.
    These are recorded IndexMap extents, not inferred nominal magnet lengths.
    """
    frame = read_ascii_sdds(path)
    if frame.attrs['units'].get('s') != 'm' or 'element_names' not in frame:
        raise ValueError('element positions require s in m and element_names')
    if len(frame) < 2 or not np.isfinite(frame.s).all() or not frame.s.is_monotonic_increasing:
        raise ValueError('element positions must be finite and nondecreasing')
    by_name = {}
    rows = list(frame[['s', 'element_names']].itertuples(index=False, name=None))
    for (start, names), (stop, _) in zip(rows[:-1], rows[1:]):
        if stop == start:
            continue
        for name in {part.strip() for part in names.split(',') if part.strip()}:
            spans = by_name.setdefault(name, [])
            if spans and spans[-1][1] == start:
                spans[-1][1] = stop
            else:
                spans.append([start, stop])
    records = [dict(element=name, s_begin_m=start, s_end_m=stop)
               for name, spans in by_name.items() for start, stop in spans]
    if not records:
        raise ValueError('no finite named element spans')
    # Longest first for coincident starts; overlapping spans occupy separate lanes.
    records.sort(key=lambda row: (row['s_begin_m'], -row['s_end_m'], row['element']))
    lane_ends = []
    for row in records:
        lane = next((i for i, end in enumerate(lane_ends) if end <= row['s_begin_m']), len(lane_ends))
        if lane == len(lane_ends):
            lane_ends.append(row['s_end_m'])
        else:
            lane_ends[lane] = row['s_end_m']
        row['lane'] = lane
    return pd.DataFrame(records)


def draw_element_strip(axis, spans):
    for row in spans.itertuples(index=False):
        y = -row.lane
        axis.hlines(y, row.s_begin_m, row.s_end_m, color='.3', lw=1.2)
        axis.vlines([row.s_begin_m, row.s_end_m], y - .12, y + .12, color='.3', lw=1.2)
        axis.text((row.s_begin_m + row.s_end_m) / 2, y + .3, row.element, ha='center', va='center',
                  fontsize=11)
    axis.set_ylim(-spans.lane.max() - .3, .7)
    axis.set_ylabel('Elements', fontsize=10)
    axis.set_yticks([])
    axis.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    for spine in axis.spines.values():
        spine.set_visible(False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).with_name('map-2-dba.in'))
    parser.add_argument('--executable', type=Path, default=Path('omp-build/src/opalx'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('field-history'))
    parser.add_argument('--plot-only', action='store_true')
    parser.add_argument('--x-coordinate', choices=['t', 's'], default='t')
    parser.add_argument('--element-positions', type=Path,
                        help='ElementPositions.sdds for an aligned element strip (requires --x-coordinate s)')
    parser.add_argument('--dt', type=float, default=1e-12)
    parser.add_argument('--maxsteps', type=int, default=20000)
    parser.add_argument('--preserve-map-settings', action='store_true',
                        help='Keep the source map options instead of selecting RK4/L1')
    args = parser.parse_args()
    if not np.isfinite(args.dt) or args.dt <= 0 or args.maxsteps < 2:
        parser.error('DT must be finite and positive; MAXSTEPS must be at least two')
    if args.element_positions is not None and args.x_coordinate != 's':
        parser.error('--element-positions requires --x-coordinate s')
    output = args.output.resolve()
    case = output / 'dba'
    case.mkdir(parents=True, exist_ok=True)
    stat = case / 'map-2-dba.stat'
    if not args.plot_only:
        if stat.exists():
            raise RuntimeError('preserve existing run: use --plot-only or another --output')
        template = args.source.read_text()
        text = template if args.preserve_map_settings else configure_map(template, .003, 1, 'RK4')
        method = re.search(r'LINEARTRANSFERMAPINTEGRATOR\s*=\s*"([^"]+)"', text)
        map_integrator = method.group(1) if method else 'BORIS'
        for key, value in [('MAXSTEPS', str(args.maxsteps)), ('DT', f'{args.dt:.16e}'), ('STATDUMPFREQ', '1'), ('PSDUMPFREQ', '0')]:
            text, count = re.subn(rf'\b{key}\s*=\s*[^,;]+', f'{key}={value}', text)
            if count != 1:
                raise RuntimeError(f'expected one {key} setting')
        (case / 'map-2-dba.in').write_text(text)
        distribution = args.source.parent.parent / 'reference-particle.txt'
        (output / 'reference-particle.txt').write_bytes(distribution.read_bytes())
        env = os.environ | {'OMP_NUM_THREADS': '1', 'OMP_PROC_BIND': 'false'}
        executable = args.executable.resolve()
        with (case / 'map-2-dba.out').open('w') as stream:
            subprocess.run([str(executable), '--info', '2', 'map-2-dba.in'], cwd=case,
                           env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
        manifest = dict(input_sha256=digest(text.encode()), template_sha256=digest(template.encode()),
                        executable_sha256=digest(executable.read_bytes()), stat_sha256=digest(stat.read_bytes()),
                        source=str(args.source.resolve()), dt_s=args.dt, maxsteps=args.maxsteps, statdumpfreq=1,
                        map_integrator=map_integrator, production_integrator='BORIS/LF2',
                        map_settings_preserved=args.preserve_map_settings, ranks=1)
        (output / 'provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    manifest = json.loads((output / 'provenance.json').read_text())
    frame = read_stat(stat)
    columns = ['t', 's', 'Bx_ref', 'By_ref', 'Bz_ref']
    frame[columns].to_csv(output / 'field-history.csv', index=False, float_format='%.12e')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.fonttype': 'none',
                         'axes.spines.top': False, 'axes.spines.right': False})
    spans = read_element_spans(args.element_positions) if args.element_positions else None
    if spans is not None:
        fig, all_axes = plt.subplots(4, 1, sharex=True, figsize=(10, 7.5), layout='constrained',
                                    gridspec_kw={'height_ratios': [.6, 1, 1.5, 1]})
        strip, axes = all_axes[0], all_axes[1:]
        draw_element_strip(strip, spans)
        boundaries = np.unique(spans[['s_begin_m', 's_end_m']].to_numpy())
        for axis in axes:
            for boundary in boundaries:
                axis.axvline(boundary, color='.65', lw=.6, ls=':', zorder=0)
        spans.to_csv(output / 'element-spans.csv', index=False, float_format='%.12e')
    else:
        fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 7), layout='constrained')
    for axis, component, color in zip(axes, ['x', 'y', 'z'], ['#0072B2', '#D55E00', '#009E73']):
        axis.plot(frame[args.x_coordinate], frame[f'B{component}_ref'], color=color, lw=1.3)
        axis.set_ylabel(rf'$B_{component}$ [T]')
        axis.grid(alpha=.25, lw=.6)
        axis.axhline(0.0, color='0.4', lw=.6, ls=':')
        if np.max(np.abs(frame[f'B{component}_ref'])) == 0.0:
            axis.set_ylim(-.01, .01)
            axis.text(.98, .76, 'zero throughout', transform=axis.transAxes, ha='right', fontsize=10)
    axes[-1].set_xlabel(r'Path length $s$ [m]' if args.x_coordinate == 's' else 'Time [ns]')
    lower, upper = frame[args.x_coordinate].min(), frame[args.x_coordinate].max()
    if spans is not None:
        lower, upper = min(lower, spans.s_begin_m.min()), max(upper, spans.s_end_m.max())
    axes[-1].set_xlim(lower, upper)
    fig.suptitle('DBA reference magnetic field from the .stat file', fontsize=14)
    fig.get_layout_engine().set(rect=(0, .07 if spans is not None else .045, 1, .95))
    fig.text(.5, .015, f'Production reference: Boris/LF2; DT = {manifest["dt_s"] * 1e12:g} ps, '
             f'statistics every step. Map selection: {manifest["map_integrator"]}.',
             ha='center', fontsize=9)
    if spans is not None:
        overlap_note = ('no overlapping spans.' if spans.lane.max() == 0
                        else 'overlapping spans occupy separate lanes.')
        fig.text(.5, .043, 'Element strip: recorded IndexMap extents; ' + overlap_note,
                 ha='center', fontsize=9)
    name = 'magnetic-field-vs-s' if args.x_coordinate == 's' else 'magnetic-field-vs-time'
    fig.savefig(output / f'{name}.png', dpi=220, facecolor='white')
    fig.savefig(output / f'{name}.svg', facecolor='white')
    plt.close(fig)
    plot_provenance = dict(stat=str(stat), stat_sha256=digest(stat.read_bytes()),
                           x_coordinate=args.x_coordinate, script_sha256=digest(Path(__file__).read_bytes()))
    if args.element_positions is not None:
        plot_provenance.update(element_positions=str(args.element_positions.resolve()),
                               element_positions_sha256=digest(args.element_positions.read_bytes()),
                               element_extent_definition='recorded IndexMap membership intervals, including overlaps')
    (output / f'{name}-provenance.json').write_text(json.dumps(plot_provenance, indent=2) + '\n')
    print(frame[columns].agg(['min', 'max']).to_string())
    print(f'{len(frame)} samples from {stat}', flush=True)


if __name__ == '__main__':
    main()
