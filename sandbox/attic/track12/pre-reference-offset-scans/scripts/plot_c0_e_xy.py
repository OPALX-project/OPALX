#!/usr/bin/env python3
"""Plot c0 primary-bunch electric-field magnitudes in the transverse x-y plane."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_INPUT = (
    Path(__file__).resolve().parent / "c0_field_compare" / "c0_field_sample_combined.csv"
)
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "c0_field_compare"


def configure_matplotlib(output_dir: Path) -> None:
    cache_root = output_dir / ".plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg", force=True)


def plot_case(data: pd.DataFrame, case: str, output_dir: Path) -> Path:
    configure_matplotlib(output_dir)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    group = data.loc[data["case"].eq(case)].copy()
    if group.empty:
        raise ValueError(f"no rows for case {case!r}")

    x_um = group["x_m"].to_numpy() * 1.0e6
    y_um = group["y_m"].to_numpy() * 1.0e6
    e_opalx = group["E_abs_opalx"].to_numpy()
    e_manufactured = group["E_abs_analytic"].to_numpy()
    ratio = group["E_abs_ratio"].to_numpy()

    e_positive = np.concatenate([e_opalx[e_opalx > 0.0], e_manufactured[e_manufactured > 0.0]])
    e_norm = LogNorm(vmin=np.nanpercentile(e_positive, 1), vmax=np.nanpercentile(e_positive, 99))
    ratio_positive = ratio[np.isfinite(ratio) & (ratio > 0.0)]
    ratio_norm = LogNorm(
        vmin=max(np.nanpercentile(ratio_positive, 1), 1.0e-3),
        vmax=np.nanpercentile(ratio_positive, 99),
    )

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0), constrained_layout=True)
    panels = [
        (axes[0], e_opalx, e_norm, r"OPALX $|E|$ [V/m]"),
        (axes[1], e_manufactured, e_norm, r"manufactured $|E|$ [V/m]"),
        (axes[2], ratio, ratio_norm, r"OPALX / manufactured $|E|$"),
    ]
    for ax, values, norm, title in panels:
        sc = ax.scatter(
            x_um,
            y_um,
            c=values,
            s=4.0,
            cmap="viridis",
            norm=norm,
            linewidths=0.0,
            rasterized=True,
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$x$ [$\mu$m]")
        ax.set_ylabel(r"$y$ [$\mu$m]")
        ax.set_title(title)
        ax.grid(True, color="0.88", lw=0.5)
        fig.colorbar(sc, ax=ax, shrink=0.88)

    fig.suptitle(f"c0 primary bunch field sample: {case.replace('_', ' ')}", fontsize=12)
    output = output_dir / f"{case}_c0_e_xy.png"
    fig.savefig(output, dpi=240)
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input_csv)
    outputs = []
    for case in data["case"].drop_duplicates():
        outputs.append(plot_case(data, str(case), args.output_dir))
    for output in outputs:
        print(f"wrote {output}")


if __name__ == "__main__":
    main()
