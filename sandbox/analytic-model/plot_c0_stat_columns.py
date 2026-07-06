#!/usr/bin/env python3
"""Plot columns 1:6 from an OPALX c0 SDDS stat file."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STAT = (
    ROOT
    / "sandbox"
    / "Drift-Experiment"
    / "one_c0_track12_edge_t0_ip6mm_400k_50ps_dt5fs_extended_bb_nobound"
    / "spacecharge_drift_withness_c0.stat"
)
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "c0_stat_columns_1_6.png"


def configure_matplotlib(output: Path) -> None:
    cache_root = output.parent / ".plot-cache"
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    cache_root.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)


def read_sdds_stat(path: Path) -> pd.DataFrame:
    text = path.read_text(encoding="utf-8")
    columns = re.findall(r"&column\s+name=([^,\n]+),", text)
    lines = [line.strip() for line in text[text.index("&data") :].splitlines() if line.strip()]
    end_index = lines.index("&end")
    payload = lines[end_index + 4 :]

    rows = []
    for line in payload:
        values = line.split()
        if len(values) == len(columns):
            rows.append(values)

    data = pd.DataFrame(rows, columns=columns)
    for column in data.columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data


def plot_columns(data: pd.DataFrame, output: Path) -> None:
    configure_matplotlib(output)
    import matplotlib.pyplot as plt

    plot_data = data.loc[:, ["t", "s", "numParticles", "charge", "energy", "rms_x"]].copy()
    plot_data["t_ps"] = plot_data["t"] * 1.0e3
    plot_data["s_mm"] = plot_data["s"] * 1.0e3
    plot_data["charge_nC"] = plot_data["charge"] * 1.0e9
    plot_data["rms_x_um"] = plot_data["rms_x"] * 1.0e6

    figure, axes = plt.subplots(5, 1, figsize=(8.0, 9.2), sharex=True, constrained_layout=True)
    series = [
        ("s_mm", "s [mm]", "tab:blue"),
        ("numParticles", "macroparticles", "tab:green"),
        ("charge_nC", "charge [nC]", "tab:red"),
        ("energy", "energy [MeV]", "tab:purple"),
        ("rms_x_um", "rms_x [um]", "tab:orange"),
    ]

    for axis, (column, ylabel, color) in zip(axes, series, strict=True):
        axis.plot(plot_data["t_ps"], plot_data[column], color=color, lw=1.6)
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=0.25)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)

    axes[-1].set_xlabel("t [ps]")
    figure.suptitle("c0 stat columns 1:6")
    figure.savefig(output, dpi=220)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stat", type=Path, default=DEFAULT_STAT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = read_sdds_stat(args.stat)
    plot_columns(data, args.output)
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
