#!/usr/bin/env python3
"""Compare the 10-step one-source OPALX BeamBeam run to the manufactured model."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TRACK12_MODULE = ROOT / "sandbox" / "track12particles" / "track12particles.py"
DEFAULT_INITIAL = ROOT / "sandbox" / "track12particles" / "outputs" / "track12_initial_conditions.csv"
DEFAULT_OPALX_DIR = ROOT / "sandbox" / "track12particles" / "opalx"
DEFAULT_PREFIX = "track12_beambeam_one_source_1fs_window6mm_binned_400k_timing_10steps"
DEFAULT_OUTPUT = DEFAULT_OPALX_DIR / "track12_one_source_10step_phase_space_comparison.csv"


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TRACK12_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_fromfile(path: Path, species: str) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as stream:
        count = int(stream.readline().strip())
        header = stream.readline().split()
        rows = [dict(zip(header, map(float, line.split()), strict=True)) for line in stream if line.strip()]
    if len(rows) != count:
        raise ValueError(f"{path}: expected {count} particles, got {len(rows)}")
    frame = pd.DataFrame(rows)
    frame["pair"] = np.arange(1, len(frame) + 1)
    frame["species"] = species
    return frame


def initial_input_deltas(initial: pd.DataFrame, opalx_dir: Path) -> pd.DataFrame:
    files = {
        "electron": opalx_dir / "track12_electrons.fromfile",
        "positron": opalx_dir / "track12_positrons.fromfile",
    }
    rows: list[dict[str, float | int | str]] = []
    for species, path in files.items():
        opalx = load_fromfile(path, species)
        ref = initial.loc[initial["species"] == species].sort_values("pair").reset_index(drop=True)
        for index, row in opalx.iterrows():
            reference = ref.iloc[index]
            rows.append(
                {
                    "species": species,
                    "pair": int(row["pair"]),
                    "dx_m": float(row["x"] - reference["x"]),
                    "dy_m": float(row["y"] - reference["y"]),
                    "ds_m": float(row["z"] - reference["s"]),
                    "dpx": float(row["px"] - reference["Px_beta_gamma"]),
                    "dpy": float(row["py"] - reference["Py_beta_gamma"]),
                    "dpz": float(row["pz"] - reference["Ps_beta_gamma"]),
                }
            )
    return pd.DataFrame(rows)


def load_opalx_final(opalx_dir: Path, prefix: str, ip_s_m: float, species_assignment: str) -> pd.DataFrame:
    if species_assignment == "h5-charge":
        species_by_container = {1: "electron", 2: "positron"}
    elif species_assignment == "cain-branch":
        species_by_container = {1: "positron", 2: "electron"}
    else:
        raise ValueError(f"unsupported species assignment {species_assignment!r}")

    rows: list[dict[str, float | int | str]] = []
    for container in (1, 2):
        path = opalx_dir / f"{prefix}_c{container}.h5"
        with h5py.File(path, "r") as h5:
            step_name = sorted(h5.keys(), key=lambda key: int(key.split("#")[1]))[-1]
            group = h5[step_name]
            ref_z = float(group.attrs["RefPartR"][2])
            time_s = float(group.attrs["TIME"][0])
            order = np.argsort(group["z"][:])
            for pair, index in enumerate(order, start=1):
                rows.append(
                    {
                        "model": "opalx",
                        "species_assignment": species_assignment,
                        "container": container,
                        "species": species_by_container[container],
                        "pair": pair,
                        "step": int(step_name.split("#")[1]),
                        "time_s": time_s,
                        "x_m": float(group["x"][index]),
                        "y_m": float(group["y"][index]),
                        "s_m": ref_z + float(group["z"][index]) - ip_s_m,
                        "px": float(group["px"][index]),
                        "py": float(group["py"][index]),
                        "pz": float(group["pz"][index]),
                    }
                )
    return pd.DataFrame(rows)


def manufactured_after_steps(
    initial: pd.DataFrame,
    module,
    source_selection: str,
    clock: str,
    steps: int,
    dt_s: float,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for _, first in initial.sort_values(["species", "pair"]).iterrows():
        position = first[["x", "y", "s"]].to_numpy(dtype=float)
        momentum = first[["Px_beta_gamma", "Py_beta_gamma", "Ps_beta_gamma"]].to_numpy(dtype=float)
        charge_units = float(first["charge_units"])
        if clock == "opalx-lab":
            t0_s = 0.0
            source_center_offset_m = 0.0
        elif clock == "cain-worldline":
            t0_s = float(first["t"]) / module.C_LIGHT
            source_center_offset_m = 0.0
        elif clock == "slide-insertion":
            t0_s = 0.0
            source_center_offset_m = -float(first["t"])
        else:
            raise ValueError(f"unsupported clock {clock!r}")

        sources = module.make_anisotropic_sources(
            module.DEFAULT_SOURCE_KINETIC_MEV,
            module.DEFAULT_SOURCE_CHARGE_C,
            module.SIGMA_XY_M,
            module.SIGMA_XY_M,
            module.SIGMA_Z_M,
            source_selection,
            0.0,
            module.DEFAULT_MC_SOURCE_PARTICLES,
            module.DEFAULT_MC_SEED,
        )
        if clock == "slide-insertion":
            shifted_sources = []
            for source in sources:
                center = source.center_t0_m.copy()
                center[2] += source_center_offset_m
                shifted_sources.append(
                    module.RigidAnisotropicGaussianSource(
                        source.name,
                        charge_C=source.charge_C,
                        sigma_lab_m=(module.SIGMA_XY_M, module.SIGMA_XY_M, module.SIGMA_Z_M),
                        beta_z=source.beta_z,
                        center_t0_m=center,
                        t0_s=0.0,
                    )
                )
            sources = tuple(shifted_sources)

        for step in range(steps):
            position, momentum, _e_field, _b_field = module.advance_with_fields(
                position,
                momentum,
                charge_units,
                t0_s + step * dt_s,
                dt_s,
                module.anisotropic_total_lab_fields,
                sources,
            )

        rows.append(
            {
                "model": "manufactured",
                "source_selection": source_selection,
                "clock": clock,
                "source_center_offset_m": source_center_offset_m,
                "species": first["species"],
                "pair": int(first["pair"]),
                "step": steps,
                "time_s": t0_s + steps * dt_s,
                "x_m": position[0],
                "y_m": position[1],
                "s_m": position[2],
                "px": momentum[0],
                "py": momentum[1],
                "pz": momentum[2],
            }
        )
    return pd.DataFrame(rows)


def compare(opalx: pd.DataFrame, manufactured: pd.DataFrame) -> pd.DataFrame:
    merged = opalx.merge(
        manufactured,
        on=["species", "pair"],
        suffixes=("_opalx", "_manufactured"),
        validate="one_to_one",
    )
    for column in ("x_m", "y_m", "s_m", "px", "py", "pz"):
        merged[f"d{column}"] = merged[f"{column}_opalx"] - merged[f"{column}_manufactured"]
    for column in ("x_m", "y_m", "s_m", "dx_m", "dy_m", "ds_m"):
        if column in merged:
            merged[column.replace("_m", "_um")] = 1.0e6 * merged[column]
    return merged


def print_summary(initial_delta: pd.DataFrame, comparison: pd.DataFrame) -> None:
    print("Initial OPALX .fromfile minus manufactured initial conditions:")
    for column in ("dx_m", "dy_m", "ds_m", "dpx", "dpy", "dpz"):
        print(f"  max |{column}| = {initial_delta[column].abs().max():.16e}")

    print("\nFinal OPALX minus manufactured at the requested step:")
    for column in ("dx_m", "dy_m", "ds_m", "dpx", "dpy", "dpz"):
        values = comparison[column]
        scale = 1.0e6 if column.endswith("_m") else 1.0
        unit = " um" if column.endswith("_m") else ""
        print(f"  max |{column}| = {values.abs().max() * scale:.9g}{unit}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-csv", type=Path, default=DEFAULT_INITIAL)
    parser.add_argument("--opalx-dir", type=Path, default=DEFAULT_OPALX_DIR)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ip-s-m", type=float, default=0.003)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--dt-s", type=float, default=1.0e-15)
    parser.add_argument(
        "--source-selection",
        choices=("copropagating", "oncoming"),
        default="copropagating",
        help="Manufactured source direction. No-copy OPALX uses the tracked, positive-pz source.",
    )
    parser.add_argument(
        "--clock",
        choices=("opalx-lab", "cain-worldline", "slide-insertion"),
        default="opalx-lab",
        help="Source time convention used by the manufactured comparison.",
    )
    parser.add_argument(
        "--species-assignment",
        choices=("h5-charge", "cain-branch"),
        default="h5-charge",
        help="How to map OPALX witness containers to species for the final-state comparison.",
    )
    args = parser.parse_args()

    module = load_track12_module()
    initial = pd.read_csv(args.initial_csv)
    initial_delta = initial_input_deltas(initial, args.opalx_dir)
    opalx = load_opalx_final(args.opalx_dir, args.prefix, args.ip_s_m, args.species_assignment)
    manufactured = manufactured_after_steps(
        initial, module, args.source_selection, args.clock, args.steps, args.dt_s
    )
    comparison = compare(opalx, manufactured)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.output_csv, index=False)
    print_summary(initial_delta, comparison)
    print(f"\nWrote {args.output_csv}")


if __name__ == "__main__":
    main()
