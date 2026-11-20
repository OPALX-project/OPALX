#!/usr/bin/env python3
"""Compare matching particle snapshots by ID; report differences, not acceptance tolerances."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def compare(old, new):
    rows = []
    for path in sorted(old.glob('*.h5')):
        with h5py.File(path) as source:
            if not any('id' in source[group] for group in source):
                continue  # Field diagnostics are serial-only, not particle snapshots.
        target = new / path.name
        if not target.exists():
            raise ValueError(f'No matching output: {target}')
        with h5py.File(path) as a, h5py.File(target) as b:
            if set(a) != set(b):
                raise ValueError(f'Different snapshot groups: {path.name}')
            for step in a:
                if 'id' not in a[step]:
                    continue
                ga, gb = a[step], b[step]
                ia, ib = np.asarray(ga['id']), np.asarray(gb['id'])
                if len(np.unique(ia)) != len(ia) or len(np.unique(ib)) != len(ib):
                    raise ValueError(f'Duplicate IDs: {path.name}/{step}')
                if not np.array_equal(np.sort(ia), np.sort(ib)):
                    raise ValueError(f'Different IDs: {path.name}/{step}')
                oa, ob = np.argsort(ia), np.argsort(ib)
                for key in ('TIME', 'GlobalTrackStep'):
                    if not np.array_equal(ga.attrs[key], gb.attrs[key]):
                        raise ValueError(f'Different {key}: {path.name}/{step}')
                for component in ('x', 'y', 'z', 'px', 'py', 'pz', 'q', 'sp'):
                    x, y = np.asarray(ga[component])[oa], np.asarray(gb[component])[ob]
                    if not np.isfinite(x).all() or not np.isfinite(y).all():
                        raise ValueError(f'Nonfinite {component}: {path.name}/{step}')
                    delta = y - x
                    rows.append(dict(file=path.name, step=step, component=component,
                                     count=len(x), max_abs=float(np.max(np.abs(delta), initial=0)),
                                     squared_error=float(np.sum(delta**2))))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('old', type=Path)
    parser.add_argument('new', type=Path)
    args = parser.parse_args()
    frame = compare(args.old, args.new)
    frame.to_csv(args.new / 'former_smoke_phase_space.csv', index=False)
    grouped = frame.groupby(['file', 'component']).agg(
        max_abs=('max_abs', 'max'), squared_error=('squared_error', 'sum'), count=('count', 'sum'))
    grouped['rms'] = np.sqrt(grouped.squared_error / grouped['count'])
    grouped.to_csv(args.new / 'former_smoke_summary.csv')
    old_manifest = json.loads((args.old / 'run_manifest.json').read_text())
    new_manifest = json.loads((args.new / 'run_manifest.json').read_text())
    report = dict(old=str(args.old.resolve()), new=str(args.new.resolve()),
                  analysis_equal=old_manifest['analysis'] == new_manifest['analysis'],
                  identity_clock_and_finiteness_checked=True,
                  acceptance='No numerical acceptance tolerance inferred; inspect phase-space differences.')
    (args.new / 'former_smoke_comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(grouped.to_string())
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
