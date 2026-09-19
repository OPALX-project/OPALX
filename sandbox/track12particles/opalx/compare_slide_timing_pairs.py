#!/usr/bin/env python3
"""Compare slide-timed pair-wise OPALX outputs with the manufactured model."""

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
DEFAULT_PAIR_DIR = ROOT / "sandbox" / "track12particles" / "opalx" / "timing_pairs"
DEFAULT_OUTPUT = DEFAULT_PAIR_DIR / "slide_timing_opalx_vs_manufactured_10steps.csv"
DEFAULT_KICK_OUTPUT = DEFAULT_PAIR_DIR / "slide_timing_first_kick_fields.csv"


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TRACK12_MODULE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_sources(module, row: pd.Series, source_selection: str):
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
    shifted = []
    for source in sources:
        center = source.center_t0_m.copy()
        center[2] -= float(row["t"])
        shifted.append(
            module.RigidAnisotropicGaussianSource(
                source.name,
                charge_C=source.charge_C,
                sigma_lab_m=(module.SIGMA_XY_M, module.SIGMA_XY_M, module.SIGMA_Z_M),
                beta_z=source.beta_z,
                center_t0_m=center,
                t0_s=0.0,
            )
        )
    return tuple(shifted)


def manufactured_after_steps(
    initial: pd.DataFrame,
    module,
    source_selection: str,
    steps: int,
    dt_s: float,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for _, first in initial.sort_values(["pair", "kind"]).iterrows():
        position = first[["x", "y", "s"]].to_numpy(dtype=float)
        momentum = first[["Px_beta_gamma", "Py_beta_gamma", "Ps_beta_gamma"]].to_numpy(dtype=float)
        charge_units = float(first["charge_units"])
        sources = make_sources(module, first, source_selection)
        for step in range(steps):
            position, momentum, _e_field, _b_field = module.advance_with_fields(
                position,
                momentum,
                charge_units,
                step * dt_s,
                dt_s,
                module.anisotropic_total_lab_fields,
                sources,
            )
        rows.append(
            {
                "model": "manufactured",
                "source_selection": source_selection,
                "clock": "slide-insertion",
                "species": first["species"],
                "pair": int(first["pair"]),
                "step": steps,
                "time_s": steps * dt_s,
                "x_m": position[0],
                "y_m": position[1],
                "s_m": position[2],
                "px": momentum[0],
                "py": momentum[1],
                "pz": momentum[2],
            }
        )
    return pd.DataFrame(rows)


def load_opalx_pairs(pair_dir: Path, ip_s_m: float, species_assignment: str) -> pd.DataFrame:
    if species_assignment == "h5-charge":
        species_by_container = {1: "electron", 2: "positron"}
    elif species_assignment == "cain-branch":
        species_by_container = {1: "positron", 2: "electron"}
    else:
        raise ValueError(f"unsupported species assignment {species_assignment!r}")

    rows: list[dict[str, float | int | str]] = []
    for pair in range(1, 7):
        stem = f"track12_pair{pair}_slide_timed_one_source_1fs_400k_10steps"
        for container in (1, 2):
            path = pair_dir / f"pair{pair}" / f"{stem}_c{container}.h5"
            with h5py.File(path, "r") as h5:
                step_name = sorted(h5.keys(), key=lambda key: int(key.split("#")[1]))[-1]
                group = h5[step_name]
                ref_z = float(group.attrs["RefPartR"][2])
                time_s = float(group.attrs["TIME"][0])
                rows.append(
                    {
                        "model": "opalx",
                        "species_assignment": species_assignment,
                        "container": container,
                        "species": species_by_container[container],
                        "pair": pair,
                        "step": int(step_name.split("#")[1]),
                        "time_s": time_s,
                        "x_m": float(group["x"][0]),
                        "y_m": float(group["y"][0]),
                        "s_m": ref_z + float(group["z"][0]) - ip_s_m,
                        "px": float(group["px"][0]),
                        "py": float(group["py"][0]),
                        "pz": float(group["pz"][0]),
                    }
                )
    return pd.DataFrame(rows)


def compare_final(opalx: pd.DataFrame, manufactured: pd.DataFrame) -> pd.DataFrame:
    merged = opalx.merge(
        manufactured,
        on=["species", "pair"],
        suffixes=("_opalx", "_manufactured"),
        validate="one_to_one",
    )
    for column in ("x_m", "y_m", "s_m", "px", "py", "pz"):
        merged[f"d{column}"] = merged[f"{column}_opalx"] - merged[f"{column}_manufactured"]
    return merged


def compare_first_kicks(
    pair_dir: Path,
    initial: pd.DataFrame,
    module,
    source_selection: str,
    species_assignment: str,
) -> pd.DataFrame:
    if species_assignment == "h5-charge":
        species_by_container = {1: "electron", 2: "positron"}
    elif species_assignment == "cain-branch":
        species_by_container = {1: "positron", 2: "electron"}
    else:
        raise ValueError(f"unsupported species assignment {species_assignment!r}")

    rows: list[dict[str, float | int | str]] = []
    for pair in range(1, 7):
        path = pair_dir / f"pair{pair}" / f"track12_pair{pair}_slide_timed_one_source_1fs_400k_10steps_kicks.csv"
        if not path.exists():
            continue
        kicks = pd.read_csv(path)
        kicks = kicks.loc[kicks["step"].eq(0)].copy()
        for _, kick in kicks.iterrows():
            species = species_by_container[int(kick["container"])]
            first = initial.loc[initial["pair"].eq(pair) & initial["species"].eq(species)].iloc[0]
            sources = make_sources(module, first, source_selection)
            position = np.array([kick["x_m"], kick["y_m"], kick["s_minus_ip_m"]], dtype=float)
            e_field, b_field = module.anisotropic_total_lab_fields(position, float(kick["time_s"]), sources)
            rows.append(
                {
                    "pair": pair,
                    "container": int(kick["container"]),
                    "species_assignment": species_assignment,
                    "species": species,
                    "time_s": float(kick["time_s"]),
                    "s_minus_ip_m": float(kick["s_minus_ip_m"]),
                    "source_selection": source_selection,
                    "Ex_opalx": float(kick["Ex_V_per_m"]),
                    "Ey_opalx": float(kick["Ey_V_per_m"]),
                    "Ez_opalx": float(kick["Ez_V_per_m"]),
                    "Bx_opalx": float(kick["Bx_T"]),
                    "By_opalx": float(kick["By_T"]),
                    "Bz_opalx": float(kick["Bz_T"]),
                    "E_abs_opalx": float(kick["E_abs_V_per_m"]),
                    "B_abs_opalx": float(kick["B_abs_T"]),
                    "Ex_manufactured": e_field[0],
                    "Ey_manufactured": e_field[1],
                    "Ez_manufactured": e_field[2],
                    "Bx_manufactured": b_field[0],
                    "By_manufactured": b_field[1],
                    "Bz_manufactured": b_field[2],
                    "E_abs_manufactured": float(np.linalg.norm(e_field)),
                    "B_abs_manufactured": float(np.linalg.norm(b_field)),
                    "dPx_opalx": float(kick["dPx"]),
                    "dPy_opalx": float(kick["dPy"]),
                    "dPz_opalx": float(kick["dPz"]),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["E_abs_ratio_opalx_over_manufactured"] = out["E_abs_opalx"] / out["E_abs_manufactured"]
        out["B_abs_ratio_opalx_over_manufactured"] = out["B_abs_opalx"] / out["B_abs_manufactured"]
    return out


def print_summary(final: pd.DataFrame, kicks: pd.DataFrame) -> None:
    print("Final OPALX minus manufactured:")
    for column in ("dx_m", "dy_m", "ds_m", "dpx", "dpy", "dpz"):
        values = final[column]
        scale = 1.0e6 if column.endswith("_m") else 1.0
        unit = " um" if column.endswith("_m") else ""
        print(f"  max |{column}| = {values.abs().max() * scale:.9g}{unit}")

    if kicks.empty:
        print("\nNo first-kick CSV rows found.")
        return
    print("\nFirst-kick field magnitude ratios, OPALX/manufactured:")
    for column in ("E_abs_ratio_opalx_over_manufactured", "B_abs_ratio_opalx_over_manufactured"):
        values = kicks[column].replace([np.inf, -np.inf], np.nan).dropna()
        print(f"  {column}: min={values.min():.6g}, median={values.median():.6g}, max={values.max():.6g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial-csv", type=Path, default=DEFAULT_INITIAL)
    parser.add_argument("--pair-dir", type=Path, default=DEFAULT_PAIR_DIR)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--kick-output-csv", type=Path, default=DEFAULT_KICK_OUTPUT)
    parser.add_argument("--ip-s-m", type=float, default=0.003)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--dt-s", type=float, default=1.0e-15)
    parser.add_argument("--source-selection", choices=("copropagating", "oncoming"), default="oncoming")
    parser.add_argument("--species-assignment", choices=("h5-charge", "cain-branch"), default="h5-charge")
    args = parser.parse_args()

    module = load_track12_module()
    initial = pd.read_csv(args.initial_csv)
    opalx = load_opalx_pairs(args.pair_dir, args.ip_s_m, args.species_assignment)
    manufactured = manufactured_after_steps(initial, module, args.source_selection, args.steps, args.dt_s)
    final = compare_final(opalx, manufactured)
    kicks = compare_first_kicks(
        args.pair_dir,
        initial,
        module,
        args.source_selection,
        args.species_assignment,
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.output_csv, index=False)
    if not kicks.empty:
        kicks.to_csv(args.kick_output_csv, index=False)

    print_summary(final, kicks)
    print(f"\nWrote {args.output_csv}")
    if not kicks.empty:
        print(f"Wrote {args.kick_output_csv}")


if __name__ == "__main__":
    main()
