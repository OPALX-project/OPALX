#!/usr/bin/env python3
"""
plot_spectrum.py - Plot decay-spectrum CSVs emitted by OPALX physics unit tests.

The CSV format (written by unit_tests/Physics/SpectrumTestSupport.h) is:

    # x_label: <axis label string>
    # columns: bin_low,bin_high,bin_center,density,count,analytic_pdf
    <bin_low>,<bin_high>,<bin_center>,<density>,<count>,<analytic_pdf>
    ...

Renders the histogram density as bars with the analytic PDF overlaid.

Usage:
    python plot_spectrum.py muon_michel_spectrum.csv --out muon_michel.png
    python plot_spectrum.py --all *.csv --outdir plots/
    python plot_spectrum.py --overlay mu_minus.csv mu_plus.csv --out muon_angular.png

Dependencies: numpy, matplotlib (no project-specific deps).
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_csv(path: Path):
    """Return (x_label, bin_low, bin_high, bin_center, density, count, analytic)."""
    x_label = ""
    rows = []
    with path.open() as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            first = row[0]
            if first.startswith("#"):
                comment = ",".join(row).lstrip("#").strip()
                if comment.startswith("x_label:"):
                    x_label = comment.split(":", 1)[1].strip()
                continue
            rows.append(row)

    arr = np.array(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 6:
        raise ValueError(
            f"{path}: expected 6 columns, got shape {arr.shape}"
        )
    return (
        x_label,
        arr[:, 0],  # bin_low
        arr[:, 1],  # bin_high
        arr[:, 2],  # bin_center
        arr[:, 3],  # density
        arr[:, 4],  # count
        arr[:, 5],  # analytic_pdf
    )


def plot_one(path: Path, out: Path, dpi: int = 150) -> None:
    x_label, lo, hi, ctr, density, count, analytic = parse_csv(path)

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    widths = hi - lo
    ax.bar(
        ctr,
        density,
        width=widths,
        align="center",
        color="#4c72b0",
        alpha=0.55,
        edgecolor="#1f3a68",
        linewidth=0.6,
        label="sampled (density)",
    )
    ax.plot(
        ctr,
        analytic,
        color="#c0392b",
        linewidth=2.0,
        label="analytic PDF",
    )

    ax.set_xlabel(x_label or "x")
    ax.set_ylabel("density")
    ax.set_title(path.stem)
    ax.set_xlim(lo[0], hi[-1])
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)

    n_total = int(count.sum())
    ax.text(
        0.99,
        0.97,
        f"N = {n_total}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="#444",
    )

    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"wrote {out}")


def plot_overlay(paths: list[Path], out: Path, dpi: int = 150) -> None:
    """Draw several CSVs on one axis: sampled density as steps, analytic as lines.

    Used to compare spectra side by side (e.g. mu- vs mu+ angular distributions),
    each sampled histogram against its own analytic shape.
    """
    palette = ["#4c72b0", "#c0392b", "#2e8b57", "#8e44ad", "#d68910", "#17a589"]

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    x_label = ""
    lo_min, hi_max = None, None
    n_total = 0
    for i, path in enumerate(paths):
        path = Path(path)
        xl, lo, hi, ctr, density, count, analytic = parse_csv(path)
        x_label = x_label or xl
        color = palette[i % len(palette)]
        # Sampled density as a step outline (edges = bin boundaries).
        edges = np.append(lo, hi[-1])
        ax.stairs(
            density,
            edges,
            color=color,
            alpha=0.55,
            linewidth=1.4,
            label=f"{path.stem} (sampled)",
        )
        # Analytic shape as a solid line in the same color.
        ax.plot(ctr, analytic, color=color, linewidth=2.0, label=f"{path.stem} (analytic)")
        lo_min = lo[0] if lo_min is None else min(lo_min, lo[0])
        hi_max = hi[-1] if hi_max is None else max(hi_max, hi[-1])
        n_total += int(count.sum())

    ax.set_xlabel(x_label or "x")
    ax.set_ylabel("density")
    ax.set_title(out.stem)
    ax.set_xlim(lo_min, hi_max)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.99,
        0.97,
        f"N = {n_total}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="#444",
    )
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"wrote {out}")


EXAMPLES = """\
examples:
  # Single CSV -> single PNG (named explicitly):
  python plot_spectrum.py muon_michel_spectrum.csv --out muon_michel.png

  # Single CSV -> single PNG (default name muon_michel_spectrum.png):
  python plot_spectrum.py muon_michel_spectrum.csv

  # Many CSVs at once, written to a target directory:
  python plot_spectrum.py --all *.csv --outdir plots/

  # Overlay several CSVs on one axis (e.g. mu- vs mu+ angular spectra):
  python plot_spectrum.py --overlay mu_minus.csv mu_plus.csv --out muon_angular.png

  # Higher-resolution output:
  python plot_spectrum.py spectrum.csv --out spectrum.png --dpi 300

CSV format expected (written by unit_tests/Physics/SpectrumTestSupport.h):
  # x_label: <axis label string>
  # columns: bin_low,bin_high,bin_center,density,count,analytic_pdf
  <bin_low>,<bin_high>,<bin_center>,<density>,<count>,<analytic_pdf>
  ...
"""


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="plot_spectrum.py",
        description=(
            "Plot decay-spectrum CSVs emitted by OPALX physics unit tests. "
            "Renders the sampled histogram density as bars with the analytic "
            "PDF overlaid as a line."
        ),
        epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "csv",
        nargs="*",
        metavar="CSV",
        help="one or more CSV files to plot",
    )
    p.add_argument(
        "--out",
        metavar="PATH",
        help="output PNG path (single-CSV mode only; default: <csv>.png next to the CSV)",
    )
    p.add_argument(
        "--outdir",
        default=".",
        metavar="DIR",
        help="output directory when plotting multiple CSVs (default: current directory)",
    )
    p.add_argument(
        "--all",
        action="store_true",
        help="treat every positional argument as a CSV and emit one PNG each into --outdir",
    )
    p.add_argument(
        "--overlay",
        action="store_true",
        help="draw all positional CSVs on a single axis (sampled steps + analytic lines)",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        metavar="N",
        help="output resolution in dots per inch (default: 150)",
    )
    args = p.parse_args(argv)

    if not args.csv:
        p.error("at least one CSV file is required")

    if args.overlay:
        out = Path(args.out) if args.out else Path(args.outdir) / "overlay.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        plot_overlay([Path(c) for c in args.csv], out, dpi=args.dpi)
        return 0

    if args.all or len(args.csv) > 1:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for csv_path in args.csv:
            csv_path = Path(csv_path)
            out = outdir / (csv_path.stem + ".png")
            plot_one(csv_path, out, dpi=args.dpi)
        return 0

    csv_path = Path(args.csv[0])
    out = Path(args.out) if args.out else csv_path.with_suffix(".png")
    plot_one(csv_path, out, dpi=args.dpi)
    return 0


if __name__ == "__main__":
    sys.exit(main())
