#!/usr/bin/env python3
"""Plot OPALX and manufactured witness field magnitudes at the tracked particles."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ELECTRON_COLOR = "#c00000"
POSITRON_COLOR = "#0017b8"
OPALX_MODEL = "OPALX"
MANUFACTURED_MODEL = "manufactured anisotropic Gaussian"


def add_elapsed_time(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "time_s" in out.columns:
        out["elapsed_time_fs"] = 1.0e15 * (
            out["time_s"] - out.groupby(["model", "species", "pair"])["time_s"].transform("min")
        )
    return out


def load_opalx_fields(path: Path, species_assignment: str, pair_mode: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {
        "step",
        "time_s",
        "container",
        "id",
        "z_m",
        "s_minus_ip_m",
        "E_abs_V_per_m",
        "B_abs_T",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")

    if "stage" in frame.columns:
        frame = frame.loc[frame["stage"] == "gather:reference-after-field"].copy()

    if species_assignment == "h5-charge":
        species_by_container = {1: "electron", 2: "positron"}
        kind_by_container = {1: 2, 2: 3}
    else:
        species_by_container = {2: "electron", 1: "positron"}
        kind_by_container = {2: 2, 1: 3}

    frame["species"] = frame["container"].map(species_by_container)
    frame["kind"] = frame["container"].map(kind_by_container)
    frame = frame.dropna(subset=["species"]).copy()

    if pair_mode == "id":
        pair_lookup = {}
        for container, group in frame.groupby("container", sort=False):
            ids = sorted(int(value) for value in group["id"].unique())
            pair_lookup.update({(container, particle_id): index for index, particle_id in enumerate(ids, 1)})
        frame["pair"] = [
            pair_lookup[(int(container), int(particle_id))]
            for container, particle_id in zip(frame["container"], frame["id"], strict=True)
        ]
    else:
        frame = frame.sort_values(["container", "step", "z_m", "id"]).copy()
        frame["pair"] = frame.groupby(["container", "step"], sort=False).cumcount() + 1

    out = frame.loc[
        :,
        [
            "step",
            "time_s",
            "container",
            "id",
            "species",
            "kind",
            "pair",
            "s_minus_ip_m",
            "E_abs_V_per_m",
            "B_abs_T",
        ],
    ].copy()
    out["model"] = OPALX_MODEL
    out["s_mm"] = 1.0e3 * out["s_minus_ip_m"]
    return add_elapsed_time(out)


def load_manufactured_fields(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"step", "time_s", "species", "kind", "pair", "s_mm", "E_abs_V_per_m", "B_abs_T"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")
    out = frame.loc[
        :,
        ["step", "time_s", "species", "kind", "pair", "s_mm", "E_abs_V_per_m", "B_abs_T"],
    ].copy()
    out["model"] = MANUFACTURED_MODEL
    return add_elapsed_time(out)


def filter_elapsed(frame: pd.DataFrame, max_elapsed_time_fs: float | None) -> pd.DataFrame:
    if max_elapsed_time_fs is None:
        return frame
    if "elapsed_time_fs" not in frame.columns:
        raise ValueError(f"{frame['model'].iloc[0]} frame has no elapsed_time_fs column")
    return frame.loc[frame["elapsed_time_fs"] <= max_elapsed_time_fs].copy()


def plot_quantity(
    frames: list[pd.DataFrame],
    quantity: str,
    output: Path,
    title: str,
    ylabel: str,
    xlim: tuple[float, float] | None,
    yscale: str,
) -> None:
    data = pd.concat(frames, ignore_index=True)
    species_colors = {"electron": ELECTRON_COLOR, "positron": POSITRON_COLOR}
    model_styles = {OPALX_MODEL: "-", MANUFACTURED_MODEL: "-."}
    model_widths = {OPALX_MODEL: 1.55, MANUFACTURED_MODEL: 1.25}

    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    for (model, species, pair), group in data.groupby(["model", "species", "pair"], sort=False):
        group = group.sort_values("s_mm")
        ax.plot(
            group["s_mm"],
            group[quantity],
            color=species_colors.get(str(species), "black"),
            ls=model_styles.get(str(model), "-"),
            lw=model_widths.get(str(model), 1.2),
            alpha=0.9,
        )

    ax.set_xlabel("S - S_IP [mm]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if yscale != "linear":
        ax.set_yscale(yscale)

    species_handles = [
        plt.Line2D([0], [0], color=ELECTRON_COLOR, lw=2.0, label="e-"),
        plt.Line2D([0], [0], color=POSITRON_COLOR, lw=2.0, label="e+"),
    ]
    model_handles = [
        plt.Line2D([0], [0], color="0.2", lw=model_widths[name], ls=model_styles[name], label=name)
        for name in model_styles
        if name in set(data["model"])
    ]
    first_legend = ax.legend(handles=species_handles, loc="upper left", fontsize=8, title="species")
    ax.add_artist(first_legend)
    ax.legend(handles=model_handles, loc="upper right", fontsize=8, title="model")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx-fields-csv", type=Path, required=True)
    parser.add_argument("--manufactured-csv", type=Path, required=True)
    parser.add_argument(
        "--e-output",
        type=Path,
        default=Path("track12_opalx_manufactured_E_abs_vs_s.png"),
    )
    parser.add_argument(
        "--b-output",
        type=Path,
        default=Path("track12_opalx_manufactured_B_abs_vs_s.png"),
    )
    parser.add_argument(
        "--opalx-species-assignment",
        choices=("cain-branch", "h5-charge"),
        default="cain-branch",
    )
    parser.add_argument("--pair-mode", choices=("z-order", "id"), default="z-order")
    parser.add_argument("--max-elapsed-time-fs", type=float)
    parser.add_argument("--xlim", nargs=2, type=float, metavar=("MIN_MM", "MAX_MM"))
    parser.add_argument("--yscale", choices=("linear", "log"), default="linear")
    args = parser.parse_args()

    opalx = filter_elapsed(
        load_opalx_fields(args.opalx_fields_csv, args.opalx_species_assignment, args.pair_mode),
        args.max_elapsed_time_fs,
    )
    manufactured = filter_elapsed(load_manufactured_fields(args.manufactured_csv), args.max_elapsed_time_fs)
    xlim = tuple(args.xlim) if args.xlim is not None else None

    plot_quantity(
        [opalx, manufactured],
        "E_abs_V_per_m",
        args.e_output,
        "|E| at witness particles",
        "|E| [V/m]",
        xlim,
        args.yscale,
    )
    plot_quantity(
        [opalx, manufactured],
        "B_abs_T",
        args.b_output,
        "|B| at witness particles",
        "|B| [T]",
        xlim,
        args.yscale,
    )
    print(args.e_output.resolve())
    print(args.b_output.resolve())
    print(
        "OPALX |E| range:",
        f"{opalx['E_abs_V_per_m'].min():.6g}",
        "to",
        f"{opalx['E_abs_V_per_m'].max():.6g}",
        "V/m",
    )
    print(
        "manufactured |E| range:",
        f"{manufactured['E_abs_V_per_m'].min():.6g}",
        "to",
        f"{manufactured['E_abs_V_per_m'].max():.6g}",
        "V/m",
    )
    print("OPALX |B| range:", f"{opalx['B_abs_T'].min():.6g}", "to", f"{opalx['B_abs_T'].max():.6g}", "T")
    print(
        "manufactured |B| range:",
        f"{manufactured['B_abs_T'].min():.6g}",
        "to",
        f"{manufactured['B_abs_T'].max():.6g}",
        "T",
    )


if __name__ == "__main__":
    main()
