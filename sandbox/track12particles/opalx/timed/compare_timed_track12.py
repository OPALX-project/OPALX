#!/usr/bin/env python3
"""Compare exact-grid OPALX witness trajectories with TestParticleOrbit.dat."""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import os
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REFERENCE = ROOT / "sandbox" / "TestParticleOrbit.dat"
TRACK12_SCRIPT = ROOT / "sandbox" / "track12particles" / "track12particles.py"
C_LIGHT = 299_792_458.0
CAIN_CT_STEP_M = 1.8e-6
IP_S_M = 4.0e-3
COLORS = {"electron": "#b51f2e", "positron": "#1f5fb5"}
LINE_STYLES = ["-", (0, (8, 6)), (0, (5, 5)), (0, (10, 8)), (0, (2, 3)), (0, (9, 4, 2, 4))]


def load_track12_module():
    spec = importlib.util.spec_from_file_location("track12particles", TRACK12_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {TRACK12_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scalar_attribute(group: h5py.Group, name: str) -> float:
    return float(np.asarray(group.attrs[name]).reshape(-1)[0])


def reference_position(group: h5py.Group) -> np.ndarray:
    if "RefPartR" in group.attrs:
        return np.asarray(group.attrs["RefPartR"], dtype=float)
    return np.array([0.0, 0.0, scalar_attribute(group, "SPOS")])


def transverse_periods(deck: Path, mesh: list[int]) -> tuple[float, float]:
    text = deck.read_text(encoding="utf-8")
    match = re.search(
        r'APERTURE\s*=\s*"RECTANGLE\(\s*([^,]+),\s*([^)]+)\)"', text
    )
    if match is None:
        raise ValueError(f"{deck}: expected a rectangular BeamBeam aperture")
    widths = (float(match.group(1)), float(match.group(2)))
    # RMin/RMax are the first and last cell centres.  The periodic particle
    # layout spans one additional cell spacing, N/(N-1) times that interval.
    return tuple(width * n / (n - 1) for width, n in zip(widths, mesh[:2], strict=True))


def periodic_delta(delta: np.ndarray, periods_m: tuple[float, float]) -> np.ndarray:
    adjusted = delta.copy()
    for axis, period in enumerate(periods_m):
        adjusted[axis] -= np.rint(adjusted[axis] / period) * period
    return adjusted


def load_opalx(
    path: Path,
    species: str,
    pair_t0_s: float,
    births: dict[int, float],
    periods_m: tuple[float, float],
) -> pd.DataFrame:
    raw_rows: list[dict[str, float | int | str]] = []
    previous: dict[int, dict[str, float | int | str]] = {}
    previous_time_s: float | None = None
    birth_order = sorted(births)
    with h5py.File(path, "r") as h5:
        groups = sorted(h5.values(), key=lambda group: int(group.name.rsplit("#", 1)[1]))
        for group in groups:
            time_s = scalar_attribute(group, "TIME")
            global_step = int(scalar_attribute(group, "GlobalTrackStep"))
            ref_r = reference_position(group)
            ct_m = C_LIGHT * (time_s - pair_t0_s)
            current: list[dict[str, float | int | str]] = []
            for index, particle_id_raw in enumerate(group["id"][:]):
                particle_id = int(particle_id_raw)
                current.append(
                    {
                        "species": species,
                        "sample_origin": "h5",
                        "particle_id": particle_id,
                        "global_step": global_step,
                        "time_s": time_s,
                        "ct_m": ct_m,
                        "x_opalx_m": ref_r[0] + float(group["x"][index]),
                        "y_opalx_m": ref_r[1] + float(group["y"][index]),
                        "x_h5_wrapped_m": ref_r[0] + float(group["x"][index]),
                        "y_h5_wrapped_m": ref_r[1] + float(group["y"][index]),
                        "s_opalx_m": ref_r[2] + float(group["z"][index]) - IP_S_M,
                        "px_opalx": float(group["px"][index]),
                        "py_opalx": float(group["py"][index]),
                        "pz_opalx": float(group["pz"][index]),
                    }
                )

            if len(current) > len(birth_order):
                raise ValueError(
                    f"{path}: step {global_step} has {len(current)} particles, "
                    f"but only {len(birth_order)} births are defined"
                )
            # The H5 sample at an exact birth time precedes insertion.  A newborn
            # therefore first appears one CAIN interval later.  The population is
            # monotonic in this no-deletion test, so its size identifies how many
            # births have occurred without relying on unstable MPI particle IDs.
            expected_pairs = birth_order[: len(current)]

            assignment: dict[int, int] = {}
            if previous:
                if previous_time_s is None:
                    raise AssertionError("previous time is missing")
                existing_pairs = sorted(previous)
                if len(current) not in (len(existing_pairs), len(existing_pairs) + 1):
                    raise ValueError(
                        f"{path}: population changed from {len(existing_pairs)} to "
                        f"{len(current)} at step {global_step}"
                    )
                dt_s = time_s - previous_time_s
                cost = np.empty((len(existing_pairs), len(current)), dtype=float)
                for row_index, pair in enumerate(existing_pairs):
                    prior = previous[pair]
                    momentum = np.array(
                        [prior["px_opalx"], prior["py_opalx"], prior["pz_opalx"]],
                        dtype=float,
                    )
                    gamma = np.sqrt(1.0 + np.dot(momentum, momentum))
                    prior_position = np.array(
                        [prior["x_opalx_m"], prior["y_opalx_m"], prior["s_opalx_m"]],
                        dtype=float,
                    )
                    predicted_position = prior_position + C_LIGHT * dt_s * momentum / gamma
                    prior_momentum = momentum
                    for column_index, candidate in enumerate(current):
                        candidate_position = np.array(
                            [
                                candidate["x_h5_wrapped_m"],
                                candidate["y_h5_wrapped_m"],
                                candidate["s_opalx_m"],
                            ],
                            dtype=float,
                        )
                        candidate_momentum = np.array(
                            [
                                candidate["px_opalx"],
                                candidate["py_opalx"],
                                candidate["pz_opalx"],
                            ],
                            dtype=float,
                        )
                        position_residual = periodic_delta(
                            candidate_position - predicted_position, periods_m
                        )
                        position_term = np.sum(np.square(position_residual / 1.0e-6))
                        momentum_term = np.sum(
                            np.square((candidate_momentum - prior_momentum) / 5.0e-2)
                        )
                        cost[row_index, column_index] = position_term + momentum_term

                best_columns: tuple[int, ...] | None = None
                best_cost = np.inf
                for columns in itertools.permutations(range(len(current)), len(existing_pairs)):
                    candidate_cost = sum(cost[row, column] for row, column in enumerate(columns))
                    if candidate_cost < best_cost:
                        best_cost = candidate_cost
                        best_columns = columns
                if best_columns is None:
                    raise AssertionError("continuity assignment has no solution")
                assignment.update(zip(existing_pairs, best_columns, strict=True))

            assigned_columns = set(assignment.values())
            new_pairs = sorted(set(expected_pairs) - set(previous))
            new_columns = sorted(set(range(len(current))) - assigned_columns)
            if len(new_pairs) != len(new_columns):
                raise ValueError(
                    f"{path}: could not identify new birth at step {global_step}: "
                    f"pairs={new_pairs}, rows={new_columns}"
                )
            assignment.update(zip(new_pairs, new_columns, strict=True))

            next_previous: dict[int, dict[str, float | int | str]] = {}
            for pair, column in sorted(assignment.items()):
                row = current[column]
                if pair in previous:
                    prior = previous[pair]
                    wrapped_delta = periodic_delta(
                        np.array(
                            [
                                row["x_h5_wrapped_m"] - prior["x_h5_wrapped_m"],
                                row["y_h5_wrapped_m"] - prior["y_h5_wrapped_m"],
                                0.0,
                            ],
                            dtype=float,
                        ),
                        periods_m,
                    )
                    row["x_opalx_m"] = prior["x_opalx_m"] + wrapped_delta[0]
                    row["y_opalx_m"] = prior["y_opalx_m"] + wrapped_delta[1]
                row["pair"] = pair
                raw_rows.append(row)
                next_previous[pair] = row
            previous = next_previous
            previous_time_s = time_s

    if set(previous) != set(births):
        raise ValueError(f"{path}: final active pairs are {sorted(previous)}, expected {sorted(births)}")
    frame = pd.DataFrame(raw_rows)
    frame["reference_step"] = [
        int(np.rint((ct_m - births[int(pair)]) / CAIN_CT_STEP_M))
        for ct_m, pair in zip(frame["ct_m"], frame["pair"], strict=True)
    ]
    residual = np.array(
        [
            ct_m - (births[int(pair)] + step * CAIN_CT_STEP_M)
            for ct_m, pair, step in zip(
                frame["ct_m"], frame["pair"], frame["reference_step"], strict=True
            )
        ]
    )
    if np.max(np.abs(residual)) > 2.0e-12:
        raise ValueError(f"{path}: OPALX times do not lie on the CAIN grid")
    return frame


def add_missing_birth_samples(
    opalx: pd.DataFrame, reference: pd.DataFrame, pair_t0_s: float
) -> pd.DataFrame:
    """Add exact input states when H5 starts one step after an exact-grid birth."""
    additions: list[dict[str, float | int | str]] = []
    present = set(
        zip(opalx["species"], opalx["pair"], opalx["reference_step"], strict=True)
    )
    for _, row in reference.loc[reference["step"] == 0].iterrows():
        key = (row["species"], int(row["pair"]), 0)
        if key in present:
            continue
        additions.append(
            {
                "species": row["species"],
                "sample_origin": "emittedfromfile input",
                "particle_id": -1,
                "global_step": -1,
                "time_s": pair_t0_s + float(row["t"]) / C_LIGHT,
                "ct_m": float(row["t"]),
                "x_opalx_m": float(row["x"]),
                "y_opalx_m": float(row["y"]),
                "x_h5_wrapped_m": float(row["x"]),
                "y_h5_wrapped_m": float(row["y"]),
                "s_opalx_m": float(row["s"]),
                "px_opalx": float(row["Px_beta_gamma"]),
                "py_opalx": float(row["Py_beta_gamma"]),
                "pz_opalx": float(row["Ps_beta_gamma"]),
                "pair": int(row["pair"]),
                "reference_step": 0,
            }
        )
    if additions:
        opalx = pd.concat([opalx, pd.DataFrame(additions)], ignore_index=True)
    return opalx


def merge_reference(reference: pd.DataFrame, opalx: pd.DataFrame) -> pd.DataFrame:
    selected = reference.rename(
        columns={
            "step": "reference_step",
            "t": "ct_cain_m",
            "x": "x_cain_m",
            "y": "y_cain_m",
            "s": "s_cain_m",
            "Px_beta_gamma": "px_cain",
            "Py_beta_gamma": "py_cain",
            "Ps_beta_gamma": "pz_cain",
        }
    )[
        [
            "species",
            "pair",
            "reference_step",
            "ct_cain_m",
            "x_cain_m",
            "y_cain_m",
            "s_cain_m",
            "px_cain",
            "py_cain",
            "pz_cain",
        ]
    ]
    merged = opalx.merge(
        selected,
        on=["species", "pair", "reference_step"],
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != len(reference):
        raise ValueError(
            f"matched {len(merged)} OPALX/reference rows, expected all {len(reference)} rows"
        )
    for coordinate in ("x", "y", "s"):
        merged[f"d{coordinate}_m"] = merged[f"{coordinate}_opalx_m"] - merged[f"{coordinate}_cain_m"]
    for coordinate in ("px", "py", "pz"):
        merged[f"d{coordinate}"] = merged[f"{coordinate}_opalx"] - merged[f"{coordinate}_cain"]
    return merged.sort_values(["species", "pair", "reference_step"]).reset_index(drop=True)


def error_metrics(comparison: pd.DataFrame) -> dict[str, object]:
    metrics: dict[str, object] = {"matched_samples": int(len(comparison))}
    for coordinate, unit_scale in (("x", 1.0e6), ("y", 1.0e6), ("s", 1.0e6)):
        difference = comparison[f"d{coordinate}_m"].to_numpy()
        reference = comparison[f"{coordinate}_cain_m"].to_numpy()
        metrics[f"rmse_{coordinate}_um"] = float(
            np.sqrt(np.mean(difference**2)) * unit_scale
        )
        metrics[f"max_abs_{coordinate}_um"] = float(
            np.max(np.abs(difference)) * unit_scale
        )
        metrics[f"relative_l2_{coordinate}"] = float(
            np.linalg.norm(difference) / np.linalg.norm(reference)
        )
    return metrics


def make_summary(
    comparison: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object], dict[str, dict[str, object]]]:
    rows: list[dict[str, float | int | str]] = []
    for (species, pair), group in comparison.groupby(["species", "pair"], sort=True):
        last = group.iloc[-1]
        rows.append(
            {
                "species": species,
                "pair": int(pair),
                "samples": len(group),
                "rmse_x_um": float(np.sqrt(np.mean(np.square(group["dx_m"]))) * 1.0e6),
                "rmse_y_um": float(np.sqrt(np.mean(np.square(group["dy_m"]))) * 1.0e6),
                "rmse_s_um": float(np.sqrt(np.mean(np.square(group["ds_m"]))) * 1.0e6),
                "final_dx_um": float(last["dx_m"] * 1.0e6),
                "final_dy_um": float(last["dy_m"] * 1.0e6),
                "final_ds_um": float(last["ds_m"] * 1.0e6),
                "final_dpx": float(last["dpx"]),
                "final_dpy": float(last["dpy"]),
                "final_dpz": float(last["dpz"]),
            }
        )
    summary = pd.DataFrame(rows)
    overall = error_metrics(comparison)
    by_species = {
        str(species): error_metrics(group)
        for species, group in comparison.groupby("species", sort=True)
    }
    return summary, overall, by_species


def make_first_kick_summary(comparison: pd.DataFrame) -> pd.DataFrame:
    first_kicks = comparison.loc[comparison["reference_step"].eq(1)].copy()
    first_kicks["px_ratio_opalx_to_cain"] = (
        first_kicks["px_opalx"] / first_kicks["px_cain"]
    )
    return first_kicks[
        [
            "species",
            "pair",
            "global_step",
            "px_opalx",
            "px_cain",
            "px_ratio_opalx_to_cain",
            "py_opalx",
            "py_cain",
        ]
    ].sort_values(["species", "pair"])


def make_identity_summary(
    opalx: pd.DataFrame, periods_m: tuple[float, float]
) -> dict[str, object]:
    h5 = opalx.loc[opalx["sample_origin"].eq("h5")]
    return {
        str(species): {
            "physical_pairs": int(group["pair"].nunique()),
            "raw_h5_ids": int(group["particle_id"].nunique()),
            "raw_h5_ids_by_pair": {
                str(int(pair)): int(pair_group["particle_id"].nunique())
                for pair, pair_group in group.groupby("pair", sort=True)
            },
            "transverse_wraps_by_pair": {
                str(int(pair)): {
                    coordinate: int(
                        np.count_nonzero(
                            np.abs(
                                np.diff(
                                    pair_group.sort_values("reference_step")[
                                        f"{coordinate}_h5_wrapped_m"
                                    ]
                                )
                            )
                            > 0.5 * periods_m[axis]
                        )
                    )
                    for axis, coordinate in enumerate(("x", "y"))
                }
                for pair, pair_group in group.groupby("pair", sort=True)
            },
        }
        for species, group in h5.groupby("species", sort=True)
    }


def configure_matplotlib(output_dir: Path) -> None:
    cache = output_dir / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    cache.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def plot_trajectories(
    comparison: pd.DataFrame, output: Path, coordinate: str, slide_limits: bool
) -> None:
    configure_matplotlib(output.parent)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 11, "axes.labelsize": 15})
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.3), dpi=180, sharex=True, sharey=True)
    for ax, model in zip(axes, ("cain", "opalx"), strict=True):
        for (species, pair), group in comparison.groupby(["species", "pair"], sort=True):
            ax.plot(
                group[f"s_{model}_m"] * 1.0e3,
                group[f"{coordinate}_{model}_m"] * 1.0e6,
                color=COLORS[species],
                linestyle=LINE_STYLES[int(pair) - 1],
                linewidth=1.0,
            )
        ax.axhline(0.0, color="black", linewidth=0.7, linestyle=(0, (2, 6)))
        ax.set_xlabel("S - IP [mm]")
        ax.grid(True, color="0.9", linewidth=0.5)
        ax.set_title("CAIN reference" if model == "cain" else "OPALX exact timed births")
        if slide_limits:
            ax.set_xlim(-0.75, 1.65)
            ax.set_ylim(-9.0, 20.0)
    axes[0].set_ylabel(rf"${coordinate}$ [$\mu$m]")
    fig.legend(
        handles=[
            Line2D([0], [0], color=COLORS["electron"], label="electron"),
            Line2D([0], [0], color=COLORS["positron"], label="positron"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.935),
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        f"12 CAIN test-particle trajectories: {coordinate}(s)", y=0.995
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--run-dir", type=Path, default=HERE)
    args = parser.parse_args()

    manifest_path = args.run_dir / "results" / "preparation_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pair_t0_s = float(manifest["pair_t0_s"])
    module = load_track12_module()
    reference = module.parse_reference_file(args.reference)
    births = {
        int(pair): float(group["t"].iloc[0])
        for pair, group in reference.sort_values("step").groupby("pair")
    }
    periods_m = transverse_periods(args.run_dir / "track12_timed.in", manifest["mesh"])
    opalx = pd.concat(
        [
            load_opalx(
                args.run_dir / "track12_timed_c1.h5",
                "electron",
                pair_t0_s,
                births,
                periods_m,
            ),
            load_opalx(
                args.run_dir / "track12_timed_c2.h5",
                "positron",
                pair_t0_s,
                births,
                periods_m,
            ),
        ],
        ignore_index=True,
    )
    opalx = add_missing_birth_samples(opalx, reference, pair_t0_s)
    comparison = merge_reference(reference, opalx)
    summary, overall, by_species = make_summary(comparison)
    first_kicks = make_first_kick_summary(comparison)
    identity = {
        "periodic_particle_layout_m": list(periods_m),
        "species": make_identity_summary(opalx, periods_m),
    }

    results_dir = args.run_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(results_dir / "track12_pointwise_comparison.csv", index=False)
    summary.to_csv(results_dir / "track12_particle_summary.csv", index=False)
    first_kicks.to_csv(results_dir / "track12_first_kick_summary.csv", index=False)
    report = {
        "overall": overall,
        "by_species": by_species,
        "h5_identity": identity,
        "first_kicks": first_kicks.to_dict(orient="records"),
        "particles": summary.to_dict(orient="records"),
    }
    (results_dir / "track12_comparison.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    plot_trajectories(
        comparison,
        results_dir / "track12_cain_vs_opalx.png",
        coordinate="x",
        slide_limits=True,
    )
    plot_trajectories(
        comparison,
        results_dir / "track12_cain_vs_opalx_wide.png",
        coordinate="x",
        slide_limits=False,
    )
    plot_trajectories(
        comparison,
        results_dir / "track12_cain_vs_opalx_y.png",
        coordinate="y",
        slide_limits=False,
    )
    print(json.dumps(overall, indent=2))
    print(json.dumps(by_species, indent=2))
    print(first_kicks.to_string(index=False))
    print(summary.to_string(index=False))
    print(f"results: {results_dir.resolve()}")


if __name__ == "__main__":
    main()
