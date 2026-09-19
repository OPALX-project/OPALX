#!/usr/bin/env python3
"""Check OPALX/CAIN/manufactured run parameters for the 12-particle comparison."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


C_LIGHT = 299_792_458.0
ELEMENTARY_CHARGE_C = 1.602_176_634e-19
WITNESS_BETA = math.sqrt(3.0) / 2.0
DEFAULT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PREFIX = Path("track12_beambeam_copied_source_1fs_window6mm_binned")
DEFAULT_INPUT = Path("track12_beambeam_copied_source_1fs_window6mm_binned.in")
DEFAULT_CAIN = DEFAULT_ROOT / "sandbox" / "track12particles" / "outputs" / "track12_trajectory.csv"
DEFAULT_MANUFACTURED = (
    DEFAULT_ROOT / "sandbox" / "track12particles" / "outputs" / "track12_manufactured_trajectory.csv"
)
DEFAULT_METADATA = Path("track12_witness_metadata.csv")


def parse_scalar_assignment(text: str, name: str) -> float | None:
    match = re.search(rf"\b{name}\s*=\s*([0-9.+\-eE]+)\s*;", text)
    return float(match.group(1)) if match else None


def parse_psdumpfreq(text: str) -> int | None:
    match = re.search(r"\bPSDUMPFREQ\s*=\s*([0-9]+)", text)
    return int(match.group(1)) if match else None


def sorted_steps(h5file: h5py.File) -> list[str]:
    return sorted(h5file.keys(), key=lambda key: int(key.split("#")[1]))


def group_time_s(group: h5py.Group) -> float | None:
    if "TIME" in group.attrs:
        return float(group.attrs["TIME"][0])
    if "T" in group.attrs:
        return float(group.attrs["T"][0])
    return None


def reference_z_m(group: h5py.Group) -> float:
    if "RefPartR" in group.attrs:
        return float(group.attrs["RefPartR"][2])
    if "SPOS" in group.attrs:
        return float(group.attrs["SPOS"][0])
    return 0.0


def summarize_h5(path: Path, ip_s_m: float) -> dict[str, float | int | str | None]:
    with h5py.File(path, "r") as h5:
        steps = sorted_steps(h5)
        first = h5[steps[0]]
        last = h5[steps[-1]]
        second = h5[steps[1]] if len(steps) > 1 else None
        first_time = group_time_s(first)
        last_time = group_time_s(last)
        second_time = group_time_s(second) if second is not None else None
        saved_dt = None if first_time is None or second_time is None else second_time - first_time

        first_s = reference_z_m(first) + first["z"][:]
        last_s = reference_z_m(last) + last["z"][:]
        return {
            "file": path.name,
            "saved_groups": len(steps),
            "first_group": steps[0],
            "last_group": steps[-1],
            "time_first_s": first_time,
            "time_last_s": last_time,
            "saved_dt_s": saved_dt,
            "charge_C": float(first.attrs["CHARGE"][0]) if "CHARGE" in first.attrs else None,
            "n_first": len(first["x"]),
            "n_last": len(last["x"]),
            "s_min_first_mm": 1.0e3 * (float(np.min(first_s)) - ip_s_m),
            "s_max_first_mm": 1.0e3 * (float(np.max(first_s)) - ip_s_m),
            "s_min_last_mm": 1.0e3 * (float(np.min(last_s)) - ip_s_m),
            "s_max_last_mm": 1.0e3 * (float(np.max(last_s)) - ip_s_m),
        }


def summarize_reference(path: Path) -> dict[str, object]:
    frame = pd.read_csv(path)
    by_species = {}
    for species, group in frame.groupby("species"):
        by_species[species] = {
            "s_min_mm": float(group["s_mm"].min()),
            "s_max_mm": float(group["s_mm"].max()),
            "x_min_um": float(group["x_um"].min()),
            "x_max_um": float(group["x_um"].max()),
            "steps_per_particle": sorted(int(v) for v in group.groupby("name").size().unique()),
        }
    first_name = frame["name"].iloc[0]
    first_particle = frame.loc[frame["name"] == first_name].sort_values("step")
    ct_steps = np.diff(first_particle["t"].to_numpy(dtype=float))
    return {
        "file": path.name,
        "rows": len(frame),
        "particles": int(frame["name"].nunique()),
        "s_min_mm": float(frame["s_mm"].min()),
        "s_max_mm": float(frame["s_mm"].max()),
        "ct_step_m_median": float(np.median(ct_steps)) if len(ct_steps) else None,
        "time_step_fs_median": float(np.median(ct_steps) / C_LIGHT * 1.0e15) if len(ct_steps) else None,
        "by_species": by_species,
    }


def print_h5_summary(summary: dict[str, object], dt_track_s: float | None) -> None:
    print(f"\n{summary['file']}")
    print(f"  charge: {summary['charge_C']:.12e} C")
    print(f"  particles: first={summary['n_first']} last={summary['n_last']}")
    print(
        "  saved groups: "
        f"{summary['saved_groups']} ({summary['first_group']}..{summary['last_group']})"
    )
    print(
        "  time: "
        f"{summary['time_first_s'] * 1.0e15:.6g}..{summary['time_last_s'] * 1.0e15:.6g} fs"
    )
    if summary["saved_dt_s"] is not None:
        print(f"  saved cadence: {summary['saved_dt_s'] * 1.0e15:.6g} fs")
    if dt_track_s is not None and summary["time_last_s"] is not None:
        print(f"  inferred integration steps: {round(summary['time_last_s'] / dt_track_s)}")
    print(
        "  S-S_IP first saved: "
        f"{summary['s_min_first_mm']:.6g}..{summary['s_max_first_mm']:.6g} mm"
    )
    print(
        "  S-S_IP last saved:  "
        f"{summary['s_min_last_mm']:.6g}..{summary['s_max_last_mm']:.6g} mm"
    )


def print_reference_summary(summary: dict[str, object]) -> None:
    print(f"\n{summary['file']}")
    print(f"  rows={summary['rows']} particles={summary['particles']}")
    print(f"  S range: {summary['s_min_mm']:.6g}..{summary['s_max_mm']:.6g} mm")
    print(f"  output ct-step median: {summary['ct_step_m_median']:.12e} m")
    print(f"  output time-step median: {summary['time_step_fs_median']:.6g} fs")
    for species, values in summary["by_species"].items():
        print(
            f"  {species}: S={values['s_min_mm']:.6g}..{values['s_max_mm']:.6g} mm, "
            f"x={values['x_min_um']:.6g}..{values['x_max_um']:.6g} um, "
            f"rows/particle={values['steps_per_particle']}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--cain-csv", type=Path, default=DEFAULT_CAIN)
    parser.add_argument("--manufactured-csv", type=Path, default=DEFAULT_MANUFACTURED)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--ip-s", type=float, default=0.003)
    parser.add_argument("--target-s-mm", type=float, default=1.65)
    args = parser.parse_args()

    input_text = args.input.read_text(encoding="utf-8")
    dt_track_s = parse_scalar_assignment(input_text, "dt_track")
    max_steps = parse_scalar_assignment(input_text, "n_steps")
    psdumpfreq = parse_psdumpfreq(input_text)

    expected_source_charge = -1.25e10 * ELEMENTARY_CHARGE_C
    expected_witness_bunch_charge = 6.0 * ELEMENTARY_CHARGE_C

    print("Deck parameters")
    print(f"  input: {args.input}")
    print(f"  DT: {dt_track_s:.12e} s" if dt_track_s else "  DT: not found")
    print(f"  PSDUMPFREQ: {psdumpfreq}" if psdumpfreq else "  PSDUMPFREQ: not found")
    print(f"  MAXSTEPS: {int(max_steps)}" if max_steps else "  MAXSTEPS: not found")
    print(f"  expected source charge: {expected_source_charge:.12e} C")
    print(f"  expected witness bunch charge magnitude: {expected_witness_bunch_charge:.12e} C")
    if dt_track_s:
        print(
            "  witness advance per step: "
            f"{WITNESS_BETA * C_LIGHT * dt_track_s * 1.0e3:.12e} mm"
        )

    for suffix in ("c0", "c1", "c2"):
        print_h5_summary(summarize_h5(args.prefix.with_name(args.prefix.name + f"_{suffix}.h5"), args.ip_s), dt_track_s)

    print_reference_summary(summarize_reference(args.cain_csv))
    print_reference_summary(summarize_reference(args.manufactured_csv))

    if args.metadata.exists() and dt_track_s:
        initial = pd.read_csv(args.metadata)
        step_mm = WITNESS_BETA * C_LIGHT * dt_track_s * 1.0e3
        print(f"\nSteps needed to reach S-S_IP = {args.target_s_mm:.6g} mm")
        for row in initial.loc[initial["kind"] == 2].itertuples():
            initial_s_mm = float(row.s_m) * 1.0e3
            needed = math.ceil((args.target_s_mm - initial_s_mm) / step_mm)
            print(f"  {row.name}: initial {initial_s_mm:.6g} mm -> {needed} steps")


if __name__ == "__main__":
    main()
