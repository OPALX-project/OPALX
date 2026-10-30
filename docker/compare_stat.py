#!/usr/bin/env python3
"""Compare an OPALX .stat file against a regression test reference.

Usage: compare_stat.py <test dir> <run dir> <test name>

Checks every `stat "<column>" last|avg <tolerance>` line of
<test dir>/<test name>.rt against <test dir>/reference/<test name>.stat.
Same logic as the smoke tests in .github/workflows/cpu-serial.yml.
Runs with several beams write one file per beam (<test name>_c0.stat, ...);
then every check is applied to each of them.
"""

import glob
import os
import sys


def parse_sdds(path):
    with open(path) as f:
        lines = f.readlines()
    cols, n_params, data_end_idx, section = [], 0, None, None
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith('&') and s != '&end':
            section = s[1:]
            if section == 'parameter':
                n_params += 1
        elif s == '&end':
            if section == 'data':
                data_end_idx = i + 1
            section = None
        elif section == 'column' and s.startswith('name='):
            cols.append(s.split('name=')[1].split(',')[0].strip().strip('"'))
    rows = [
        [v.strip() for v in l.strip().split('\t')]
        for l in lines[data_end_idx + n_params:]
        if l.strip()
    ]
    return cols, rows


def stat_value(cols, rows, col, agg):
    j = cols.index(col)
    vals = [float(row[j]) for row in rows]
    if agg == 'last':
        return vals[-1]
    if agg == 'avg':
        return sum(vals) / len(vals)
    raise ValueError(f'Unknown aggregation type: {agg}')


def parse_rt(path):
    checks = []
    with open(path) as f:
        for line in f:
            line = line.split('#')[0].strip()
            parts = line.split()
            if len(parts) >= 4 and parts[0] == 'stat':
                checks.append((parts[1].strip('"'), parts[2], float(parts[3])))
    return checks


def reference_stat_files(test_dir, test):
    if os.path.exists(f'{test_dir}/reference/{test}.stat'):
        return [f'{test}.stat']
    pattern = f'{glob.escape(test_dir)}/reference/{glob.escape(test)}_c*.stat'
    return sorted(os.path.basename(p) for p in glob.glob(pattern))


def main(test_dir, run_dir, test):
    checks = parse_rt(f'{test_dir}/{test}.rt')
    stat_files = reference_stat_files(test_dir, test)
    if not stat_files:
        print(f'No reference .stat file in {test_dir}/reference')
        return 1

    failed = False
    for stat_file in stat_files:
        if not os.path.exists(f'{run_dir}/{stat_file}'):
            print(f'  [FAIL] {stat_file} was not written')
            failed = True
            continue
        ref_cols, ref_rows = parse_sdds(f'{test_dir}/reference/{stat_file}')
        got_cols, got_rows = parse_sdds(f'{run_dir}/{stat_file}')

        print(f'Checking {len(checks)} assertion(s) from {test}.rt against {stat_file} ({len(got_rows)} time steps)')
        for col, agg, tol in checks:
            r = stat_value(ref_cols, ref_rows, col, agg)
            g = stat_value(got_cols, got_rows, col, agg)
            rel_err = abs(g - r) / (abs(r) + 1e-30)
            status = 'PASS' if rel_err < tol else 'FAIL'
            if status == 'FAIL':
                failed = True
            print(f'  [{status}] {col} ({agg}): ref={r:.15e}  got={g:.15e}  rel_err={rel_err:.2e}  tol={tol:.0e}')
    return 1 if failed else 0


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
