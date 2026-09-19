#!/usr/bin/env python3
"""Plot the digitized data extracted from the provided graph.

This script can also overlay OPALX data from an ASCII SDDS/stat file.

Examples
--------
Recreate the digitized plot:
    python plot_extracted_graph.py

Overlay an OPALX stat file:
    python plot_extracted_graph.py --input-opalx drift-1.stat \
        --output replotted_graph_with_opalx.png

If your OPALX time column is stored in seconds, convert it to ns explicitly:
    python plot_extracted_graph.py --input-opalx drift-1.stat \
        --opalx-time-scale 1e9 --opalx-eps-scale 1e6

The extracted CSV is expected to contain these columns:
    t_ns, OPAL_x_m, ImpactT_x_m, OPAL_epsx_mm_mr, ImpactT_epsx_mm_mr

For OPALX stat files, the defaults are the user-requested columns:
    1:6  -> OPALX x
    1:12 -> OPALX epsilon_x
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import pandas as pd


EXPECTED_COLUMNS = [
    "t_ns",
    "OPAL_x_m",
    "ImpactT_x_m",
    "OPAL_epsx_mm_mr",
    "ImpactT_epsx_mm_mr",
]


@dataclass(frozen=True)
class SddsColumn:
    """Minimal SDDS column description."""

    name: str
    units: str | None = None
    description: str | None = None


def build_parser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    default_input = script_dir / "extracted_graph_data.csv"
    default_output = script_dir / "replotted_graph.png"

    parser = argparse.ArgumentParser(description="Replot extracted graph data.")
    parser.add_argument(
        "--input",
        type=Path,
        default=default_input,
        help="Path to the CSV file with digitized data.",
    )
    parser.add_argument(
        "--input-opalx",
        type=Path,
        default=None,
        help="Optional OPALX ASCII SDDS/stat file to overlay.",
    )
    parser.add_argument(
        "--opalx-time-col",
        type=int,
        default=1,
        help="1-based OPALX time column index. Default: 1",
    )
    parser.add_argument(
        "--opalx-x-col",
        type=int,
        default=6,
        help="1-based OPALX rms beam-size x column index. Default: 6",
    )
    parser.add_argument(
        "--opalx-eps-col",
        type=int,
        default=12,
        help="1-based OPALX epsilon_x column index. Default: 12",
    )
    parser.add_argument(
        "--opalx-time-scale",
        type=float,
        default=None,
        help=(
            "Scale factor applied to the OPALX time column. "
            "If omitted, the script auto-detects from SDDS units "
            "(1e9 for 's', 1 for 'ns')."
        ),
    )
    parser.add_argument(
        "--opalx-eps-scale",
        type=float,
        default=1.0e6,
        help="Scale factor applied to the OPALX epsilon_x column. Default: 1e6",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output,
        help="Path to save the rendered plot.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively in addition to saving it.",
    )
    return parser


def validate_columns(df: pd.DataFrame) -> None:
    missing = [col for col in EXPECTED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            "Input CSV is missing required columns: " + ", ".join(missing)
        )


def _parse_sdds_block(lines: Sequence[str], start_index: int) -> tuple[dict[str, str], int]:
    """Parse an SDDS &column or &parameter block."""
    metadata: dict[str, str] = {}
    idx = start_index + 1
    while idx < len(lines):
        line = lines[idx].strip()
        if line == "&end":
            return metadata, idx
        line = line.rstrip(",")
        if "=" in line:
            key, value = line.split("=", 1)
            metadata[key.strip()] = value.strip().strip('"')
        idx += 1
    raise ValueError("Malformed SDDS header: missing '&end'.")


def read_opalx_stat(
    path: Path,
    time_col: int = 1,
    x_col: int = 6,
    eps_col: int = 12,
    time_scale: float | None = None,
    eps_scale: float = 1.0e6,
) -> tuple[pd.DataFrame, list[SddsColumn], float, float]:
    """Read an ASCII SDDS/stat file and return transformed plotting columns.

    The parser is intentionally lightweight and only supports the ASCII layout
    used by OPALX stat files: it reads column definitions from the SDDS header
    and then extracts all numeric rows with the expected column count.
    """
    if not path.exists():
        raise FileNotFoundError(f"OPALX stat file not found: {path}")

    lines = path.read_text(errors="replace").splitlines()

    columns: list[SddsColumn] = []
    data_start: int | None = None
    idx = 0
    while idx < len(lines):
        tag = lines[idx].strip()
        if tag == "&column":
            block, idx = _parse_sdds_block(lines, idx)
            columns.append(
                SddsColumn(
                    name=block.get("name", f"col{len(columns) + 1}"),
                    units=block.get("units"),
                    description=block.get("description"),
                )
            )
        elif tag == "&parameter":
            _, idx = _parse_sdds_block(lines, idx)
        elif tag == "&data":
            _, idx = _parse_sdds_block(lines, idx)
            data_start = idx + 1
            break
        idx += 1

    if data_start is None:
        raise ValueError("Could not find SDDS '&data' section in OPALX stat file.")
    if not columns:
        raise ValueError("Could not find any SDDS column definitions in OPALX stat file.")

    ncols = len(columns)
    rows: list[list[float]] = []
    for raw_line in lines[data_start:]:
        parts = raw_line.strip().split()
        if len(parts) != ncols:
            continue
        try:
            rows.append([float(value) for value in parts])
        except ValueError:
            # Parameter/string lines immediately after &data are skipped here.
            continue

    if not rows:
        raise ValueError("No numeric data rows found in OPALX stat file.")

    raw_df = pd.DataFrame(rows, columns=[col.name for col in columns])

    if not (1 <= time_col <= ncols and 1 <= x_col <= ncols and 1 <= eps_col <= ncols):
        raise ValueError(
            f"Requested OPALX columns must be within 1..{ncols}; got "
            f"time={time_col}, x={x_col}, eps={eps_col}."
        )

    time_meta = columns[time_col - 1]
    eps_meta = columns[eps_col - 1]

    if time_scale is None:
        units = (time_meta.units or "").strip().lower()
        if units in {"s", "sec", "second", "seconds"}:
            inferred_time_scale = 1.0e9
        else:
            inferred_time_scale = 1.0
    else:
        inferred_time_scale = time_scale

    transformed = pd.DataFrame(
        {
            "t_ns": raw_df.iloc[:, time_col - 1] * inferred_time_scale,
            "OPALX_x_m": raw_df.iloc[:, x_col - 1],
            "OPALX_epsx_mm_mr": raw_df.iloc[:, eps_col - 1] * eps_scale,
        }
    )

    return transformed, columns, inferred_time_scale, eps_scale


def plot_data(
    digitized_df: pd.DataFrame | None,
    opalx_df: pd.DataFrame | None,
    output_path: Path,
    show: bool = False,
) -> None:
    if digitized_df is None and opalx_df is None:
        raise ValueError("No data available to plot.")

    fig, ax1 = plt.subplots(figsize=(8.6, 5.6), dpi=150)
    ax2 = ax1.twinx()

    if digitized_df is not None:
        ax1.plot(
            digitized_df["t_ns"],
            digitized_df["OPAL_x_m"],
            linestyle="--",
            marker="x",
            linewidth=1.0,
            markersize=5,
            color="red",
            alpha=0.7,
            label="OPAL x",
        )
        ax1.plot(
            digitized_df["t_ns"],
            digitized_df["ImpactT_x_m"],
            linestyle="--",
            marker="x",
            linewidth=1.0,
            markersize=5,
            color="limegreen",
            alpha=0.85,
            label="ImpactT x",
        )
        ax2.plot(
            digitized_df["t_ns"],
            digitized_df["OPAL_epsx_mm_mr"],
            linestyle="None",
            marker="s",
            markersize=5,
            color="red",
            label=r"OPAL $\epsilon_x$",
        )
        ax2.plot(
            digitized_df["t_ns"],
            digitized_df["ImpactT_epsx_mm_mr"],
            linestyle="None",
            marker="s",
            markersize=5,
            color="lime",
            label=r"ImpactT $\epsilon_x$",
        )

    if opalx_df is not None:
        ax1.plot(
            opalx_df["t_ns"],
            opalx_df["OPALX_x_m"],
            linestyle="-",
            linewidth=1.4,
            color="tab:blue",
            label="OPALX x",
        )
        ax2.plot(
            opalx_df["t_ns"],
            opalx_df["OPALX_epsx_mm_mr"],
            linestyle="-",
            linewidth=1.4,
            color="tab:purple",
            label=r"OPALX $\epsilon_x$",
        )

    ax1.set_xlabel("t (nS)", fontsize=14)
    ax1.set_ylabel("rms beam size [m]", fontsize=14)
    ax2.set_ylabel("rms emittance [mm-mr]", fontsize=14)

    t_max = 80.0
    left_y_max = 0.035
    right_y_max = 2.5

    if digitized_df is not None:
        t_max = max(t_max, float(digitized_df["t_ns"].max()) * 1.02)
        left_y_max = max(left_y_max, float(digitized_df[["OPAL_x_m", "ImpactT_x_m"]].max().max()) * 1.02)
        right_y_max = max(
            right_y_max,
            float(
                digitized_df[["OPAL_epsx_mm_mr", "ImpactT_epsx_mm_mr"]].max().max()
            )
            * 1.02,
        )

    if opalx_df is not None:
        t_max = max(t_max, float(opalx_df["t_ns"].max()) * 1.02)
        left_y_max = max(left_y_max, float(opalx_df["OPALX_x_m"].max()) * 1.02)
        right_y_max = max(right_y_max, float(opalx_df["OPALX_epsx_mm_mr"].max()) * 1.02)

    ax1.set_xlim(0, t_max)
    ax1.set_ylim(0, left_y_max)
    ax2.set_ylim(0, right_y_max)

    ax1.grid(True, linestyle=":", linewidth=0.6, alpha=0.75)

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(
        handles1 + handles2,
        labels1 + labels2,
        loc="lower right",
        frameon=False,
        fontsize=10,
    )

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)


def load_digitized_csv(input_path: Path, default_input: Path, allow_missing_default: bool) -> pd.DataFrame | None:
    if input_path.exists():
        df = pd.read_csv(input_path)
        validate_columns(df)
        return df

    if allow_missing_default and input_path.resolve() == default_input.resolve():
        print(
            f"Warning: default digitized CSV not found: {input_path}. Continuing without it.",
            file=sys.stderr,
        )
        return None

    raise FileNotFoundError(f"Input CSV not found: {input_path}")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    default_input = script_dir / "extracted_graph_data.csv"

    try:
        digitized_df = load_digitized_csv(
            input_path=args.input,
            default_input=default_input,
            allow_missing_default=args.input_opalx is not None,
        )

        opalx_df: pd.DataFrame | None = None
        if args.input_opalx is not None:
            opalx_df, opalx_columns, used_time_scale, used_eps_scale = read_opalx_stat(
                path=args.input_opalx,
                time_col=args.opalx_time_col,
                x_col=args.opalx_x_col,
                eps_col=args.opalx_eps_col,
                time_scale=args.opalx_time_scale,
                eps_scale=args.opalx_eps_scale,
            )
            print(
                "Loaded OPALX stat file: "
                f"time=col {args.opalx_time_col} ({opalx_columns[args.opalx_time_col - 1].name}, "
                f"units={opalx_columns[args.opalx_time_col - 1].units or 'unknown'}) x {used_time_scale:g}, "
                f"x=col {args.opalx_x_col} ({opalx_columns[args.opalx_x_col - 1].name}), "
                f"eps=col {args.opalx_eps_col} ({opalx_columns[args.opalx_eps_col - 1].name}, "
                f"units={opalx_columns[args.opalx_eps_col - 1].units or 'unknown'}) x {used_eps_scale:g}"
            )

        plot_data(digitized_df=digitized_df, opalx_df=opalx_df, output_path=args.output, show=args.show)
    except Exception as exc:  # pragma: no cover
        print(f"Error while plotting: {exc}", file=sys.stderr)
        return 1

    print(f"Saved plot to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
