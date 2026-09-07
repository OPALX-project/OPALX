#!/usr/bin/env python3
"""Portable full-map regression checker. Python standard library only.

Every entry must satisfy |M-Mref| <= atol + rtol*|Mref| after positions are scaled
by 1 metre. A nonzero absolute tolerance is essential for expected zero entries.
Missing, duplicate or nonfinite maps fail. Frozen references are never updated.
"""
import argparse
import json
import math
from pathlib import Path
import re


def parse_map(text):
    blocks = re.findall(r'Combined (?:one-turn )?linear transfer map[^\n]*\n((?:[^\n]*\n){6})', text)
    if len(blocks) != 1:
        raise ValueError(f'Expected exactly one combined map, found {len(blocks)}')
    rows = []
    for line in blocks[0].splitlines():
        # Parse every token after the logger prefix, so NaN/Inf cannot be skipped.
        values = line.split('>', 1)[-1].strip().split()
        if len(values) != 6:
            raise ValueError('Expected six entries in each matrix row')
        row = [float(v) for v in values]
        if not all(math.isfinite(v) for v in row):
            raise ValueError('Nonfinite map entry')
        rows.append(row)
    return rows


def compare(matrix, reference):
    target = reference['matrix']
    if len(matrix) != 6 or any(len(row) != 6 for row in matrix):
        raise ValueError('Measured matrix must be 6x6')
    if len(target) != 6 or any(len(row) != 6 for row in target):
        raise ValueError('Reference must be 6x6')
    atol, rtol = reference['atol'], reference['rtol']
    if not math.isfinite(atol) or not math.isfinite(rtol) or atol <= 0 or rtol < 0:
        raise ValueError('Invalid regression tolerances')
    differences, ratios = [], []
    for i in range(6):
        row, scaled = [], []
        for j in range(6):
            if not math.isfinite(target[i][j]) or not math.isfinite(matrix[i][j]):
                raise ValueError('Nonfinite reference or measured entry')
            delta = matrix[i][j]-target[i][j]
            row.append(delta)
            scaled.append(abs(delta)/(atol+rtol*abs(target[i][j])))
        differences.append(row)
        ratios.append(scaled)
    maximum = max(abs(v) for row in differences for v in row)
    ratio = max(v for row in ratios for v in row)
    worst = max(((i, j) for i in range(6) for j in range(6)), key=lambda ij: ratios[ij[0]][ij[1]])
    return dict(status='PASS' if ratio <= 1 else 'FAIL', maximum_scaled_entry_error=maximum,
                maximum_tolerance_ratio=ratio, worst_tolerance_entry=f'R{worst[0]+1}{worst[1]+1}',
                R16=matrix[0][5], R26=matrix[1][5], delta_R16=differences[0][5],
                delta_R26=differences[1][5], differences=differences, matrix=matrix)


def check_files(stdout, expected):
    return compare(parse_map(Path(stdout).read_text()), json.loads(Path(expected).read_text()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stdout', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    try:
        result = check_files(args.stdout, args.reference)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, f'FAIL: {error}\n')
    if args.json:
        args.json.write_text(json.dumps(result, indent=2)+'\n')
    print(f"{result['status']}: max error={result['maximum_scaled_entry_error']:.6e}, "
          f"tolerance ratio={result['maximum_tolerance_ratio']:.6e}, "
          f"worst={result['worst_tolerance_entry']}; all 36 entries checked")
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
