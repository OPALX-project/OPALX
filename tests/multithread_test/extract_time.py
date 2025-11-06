#!/usr/bin/env python3
import sys
import argparse
import re
import json

# This code is written by chatgpt o4 is held together by the shear hope that it will not break
def parse_wall_avg(lines):
    """
    Extract a dict { timer_name : wall_avg } from the “Wall avg” table.
    Accepts either a list of lines or one big string.
    """
    # allow passing in a big string
    if isinstance(lines, str):
        lines = lines.splitlines()

    results   = {}
    in_table  = False

    for L in lines:
        # once we see the header line, start parsing
        if not in_table:
            if 'Wall max' in L and 'Wall avg' in L:
                in_table = True
            continue

        # skip blank lines
        if not L.strip():
            continue

        # split the line into tokens
        parts = L.split()
        # we expect at least: [ name... , ranks, max, min, avg ]
        if len(parts) < 5:
            continue

        # the 4th-from-the-end token *must* be the ranks (an integer)
        if not parts[-4].isdigit():
            continue

        # grab them off the back
        avg_str  = parts[-1]
        min_str  = parts[-2]
        max_str  = parts[-3]
        ranks_str= parts[-4]

        # the rest is the name
        name = ' '.join(parts[:-4]).rstrip('.')

        # convert avg to float
        try:
            avg = float(avg_str)
        except ValueError:
            continue

        results[name] = avg

    return results
    