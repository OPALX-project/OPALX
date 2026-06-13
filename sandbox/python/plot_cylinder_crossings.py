#!/usr/bin/env python3
"""Plot first cylinder-edge crossing histograms from OPALX H5 trajectory files."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TRACK_DIR = ROOT / "sandbox/track-e-p"
DEFAULT_OUTPUT = ROOT / "sandbox/note/figs/gamma_gamma_large_cylinder_staged_dt_crossings.png"


def configure_plot_environment(output: Path) -> None:
    cache_dir = output.parent / ".plot-cache"
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))
    output.parent.mkdir(parents=True, exist_ok=True)


def first_crossing_positions(stem: str, container: str, aperture_radius_m: float) -> dict[int, float]:
    path = TRACK_DIR / f"{stem}_{container}.h5"
    first_crossing_z: dict[int, float] = {}
    with h5py.File(path, "r") as h5:
        steps = sorted(
            (key for key in h5.keys() if key.startswith("Step#")),
            key=lambda key: int(key.split("#", 1)[1]),
        )
        for step in steps:
            group = h5[step]
            ids = np.asarray(group["id"], dtype=int)
            x = np.asarray(group["x"], dtype=float)
            y = np.asarray(group["y"], dtype=float)
            z = np.asarray(group["z"], dtype=float)
            outside = np.hypot(x, y) >= aperture_radius_m
            for particle_id, z_position in zip(ids[outside], z[outside], strict=True):
                first_crossing_z.setdefault(int(particle_id), float(z_position))
    return first_crossing_z


def plot_histograms(args: argparse.Namespace) -> None:
    configure_plot_environment(args.output)
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 220, "font.size": 10})
    edges = np.linspace(args.hist_min_m, args.hist_max_m, args.hist_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = 0.58 * (edges[1] - edges[0])

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.8), sharex=True, constrained_layout=True)
    summary_lines: list[str] = []
    species = [("c1", r"$e^-$", "#2AA6B8"), ("c2", r"$e^+$", "#D33682")]

    for ax, (container, species_label, color) in zip(axes, species, strict=True):
        nominal = first_crossing_positions(args.stem, container, args.aperture_radius_m)
        nominal_values = np.asarray(list(nominal.values()), dtype=float)
        nominal_counts, _ = np.histogram(nominal_values, bins=edges)
        summary = f"{species_label}: N={int(nominal_counts.sum())}"

        if args.compare_stem:
            compare = first_crossing_positions(args.compare_stem, container, args.aperture_radius_m)
            compare_values = np.asarray(list(compare.values()), dtype=float)
            compare_counts, _ = np.histogram(compare_values, bins=edges)
            diff = compare_counts - nominal_counts
            shared_ids = sorted(set(nominal).intersection(compare))
            max_shift_um = 0.0
            if shared_ids:
                max_shift_um = 1.0e6 * max(
                    abs(compare[particle_id] - nominal[particle_id]) for particle_id in shared_ids
                )
            summary = (
                f"{species_label}: N={int(nominal_counts.sum())}/{int(compare_counts.sum())}, "
                f"|bin diff|={int(np.abs(diff).sum())}, max |dz|={max_shift_um:.1f} um"
            )
            ax.step(centers, nominal_counts, where="mid", color=color, linewidth=1.8, label=args.label)
            ax.step(
                centers,
                compare_counts,
                where="mid",
                color="black",
                linewidth=1.2,
                linestyle="--",
                label=args.compare_label,
            )
            ax.bar(
                centers,
                diff,
                width=width,
                color="0.72",
                edgecolor="0.35",
                linewidth=0.25,
                alpha=0.75,
                label="compare - nominal",
            )
            ax.axhline(0.0, color="0.25", linewidth=0.6)
        else:
            ax.bar(
                centers,
                nominal_counts,
                width=width,
                color=color,
                edgecolor="0.25",
                linewidth=0.35,
                alpha=0.82,
                label=args.label,
            )

        summary_lines.append(summary)
        ax.axvspan(args.cylinder_begin_m, args.cylinder_end_m, color="0.94", zorder=-2)
        ax.axvline(0.5 * (args.cylinder_begin_m + args.cylinder_end_m), color="0.25", linewidth=0.7)
        ax.set_ylabel(f"{species_label}\ncounts")
        ax.grid(axis="y", color="0.88", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].legend(loc="upper right", frameon=False, ncols=3, fontsize=8)
    axes[1].set_xlabel("first cylinder-edge crossing z [m]")
    fig.suptitle(args.title)
    axes[1].text(
        0.01,
        -0.34,
        "\n".join(summary_lines),
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.output}")
    for line in summary_lines:
        print(line)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", default="gamma_gamma_pairs-large-cylinder-retire1000ps-staged-dt")
    parser.add_argument("--compare-stem")
    parser.add_argument("--label", default="staged 1 fs overlap")
    parser.add_argument("--compare-label", default="comparison")
    parser.add_argument("--title", default="Staged-DT 1000 ps retirement: cylinder crossing histogram")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--aperture-radius-m", type=float, default=0.15)
    parser.add_argument("--cylinder-begin-m", type=float, default=0.01)
    parser.add_argument("--cylinder-end-m", type=float, default=0.33)
    parser.add_argument("--hist-min-m", type=float, default=-0.06)
    parser.add_argument("--hist-max-m", type=float, default=0.47)
    parser.add_argument("--hist-bins", type=int, default=18)
    return parser.parse_args()


def main() -> int:
    plot_histograms(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
