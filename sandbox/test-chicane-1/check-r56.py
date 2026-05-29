#!/usr/bin/env python3
"""Extract the R56 probe from the test-chicane-1 exit monitor H5 file."""

from __future__ import annotations

import argparse
import re

import h5py
import numpy as np


EXPECTED_R56 = -4.333333333333334e-2
INITIAL_DELTA_BY_ID = np.array([-2.0e-3, -1.0e-3, 0.0, 1.0e-3, 2.0e-3])


def step_index(name: str) -> int:
    match = re.fullmatch(r"Step#(\d+)", name)
    if not match:
        return -1
    return int(match.group(1))


def h5_attr_scalar(group: h5py.Group, name: str) -> float:
    if name not in group.attrs:
        return float("nan")
    value = np.asarray(group.attrs[name], dtype=float)
    if value.size == 0:
        return float("nan")
    return float(value.reshape(-1)[0])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5file", nargs="?", default="test-chicane-1_exit.h5")
    args = parser.parse_args()

    rows = []
    with h5py.File(args.h5file, "r") as h5:
        for name in sorted(h5.keys(), key=step_index):
            group = h5[name]
            if not {"z", "pz", "id"}.issubset(group.keys()):
                continue
            z = np.asarray(group["z"], dtype=float)
            pz = np.asarray(group["pz"], dtype=float)
            particle_id = np.asarray(group["id"], dtype=int)
            if np.any(particle_id < 0) or np.any(particle_id >= INITIAL_DELTA_BY_ID.size):
                raise RuntimeError(
                    f"{name} contains particle ids outside the configured delta scan: "
                    f"{particle_id}"
                )
            delta = INITIAL_DELTA_BY_ID[particle_id]
            slope, intercept = np.polyfit(delta, z, 1)
            spos = h5_attr_scalar(group, "SPOS")
            rows.append((name, spos, slope, intercept, particle_id, delta, z, pz))

    if not rows:
        raise RuntimeError(f"no usable Step# groups found in {args.h5file}")

    # The temporal monitor can produce an initial bookkeeping passage.  For the
    # R56 probe use the near-exit passage whose fitted reference-particle offset
    # is closest to zero.
    exit_rows = [row for row in rows if row[1] > 0.5]
    if not exit_rows:
        exit_rows = rows
    selected = min(exit_rows, key=lambda row: abs(row[3]))
    name, spos, slope, intercept, particle_id, delta, z, pz = selected
    r56_transport = -slope

    print(f"file: {args.h5file}")
    print(f"selected step: {name}")
    if np.isfinite(spos):
        print(f"s position: {spos:.12e} m")
    else:
        print("s position: not available in H5")
    print(f"fit dz_monitor/ddelta: {slope:.12e} m")
    print(f"fit intercept: {intercept:.12e} m")
    print(f"transport R56 = -dz_monitor/ddelta: {r56_transport:.12e} m")
    print(f"small-angle expected R56: {EXPECTED_R56:.12e} m")
    print(f"difference: {r56_transport - EXPECTED_R56:.12e} m")
    print()
    print("id delta z_monitor[m] pz_exit[beta*gamma]")
    for pid, d, zz, ppz in sorted(zip(particle_id, delta, z, pz)):
        print(f"{pid:d} {d:+.12e} {zz:+.12e} {ppz:.12e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
