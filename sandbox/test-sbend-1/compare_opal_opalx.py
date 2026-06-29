#!/usr/bin/env python3
"""Compare OPALX and old OPAL SBEND output for test-sbend-1.

The two executables both write a ``data`` directory, so this script expects the
two cases to be run in separate directories.  It compares only
``*_DesignPath.dat`` and ``*.stat``; phase-space HDF5 and load-balance files are
intentionally ignored.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/opalx-cache")

import numpy as np
import pandas as pd


DESIGN_COLUMNS = [
    "s",
    "Rx",
    "Ry",
    "Rz",
    "Px",
    "Py",
    "Pz",
    "Efx",
    "Efy",
    "Efz",
    "Bfx",
    "Bfy",
    "Bfz",
    "Ekin",
    "t",
]


def read_design_path(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < len(DESIGN_COLUMNS):
            continue
        rows.append([float(value) for value in fields[: len(DESIGN_COLUMNS)]])
    return pd.DataFrame(rows, columns=DESIGN_COLUMNS)


def read_sdds_stat(path: Path) -> pd.DataFrame:
    lines = path.read_text().splitlines()
    names: list[str] = []
    in_column = False
    data_start = None

    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("&column"):
            in_column = True
            continue
        if in_column and stripped.startswith("name="):
            names.append(stripped.split("=", 1)[1].split(",", 1)[0].strip())
            continue
        if in_column and stripped == "&end":
            in_column = False
            continue
        if stripped.startswith("&data"):
            for data_index in range(index + 1, len(lines)):
                if lines[data_index].strip() == "&end":
                    data_start = data_index + 1
                    break
            break

    if data_start is None:
        raise ValueError(f"Could not find SDDS data section in {path}")

    rows = []
    for line in lines[data_start:]:
        fields = line.split()
        if len(fields) < len(names):
            continue
        try:
            rows.append([float(value) for value in fields[: len(names)]])
        except ValueError:
            continue

    return pd.DataFrame(rows, columns=names)


def add_normalized_coordinates(
    opalx: pd.DataFrame, opal: pd.DataFrame, *, s_offset: float, z_offset: float
) -> tuple[pd.DataFrame, pd.DataFrame]:
    opalx = opalx.copy()
    opal = opal.copy()

    opalx["s_rel"] = opalx["s"]
    opal["s_rel"] = opal["s"] - s_offset

    if "Rz" in opalx:
        opalx["Rz_rel"] = opalx["Rz"]
        opal["Rz_rel"] = opal["Rz"] - z_offset

    if "mean_s" in opalx:
        opalx["mean_s_rel"] = opalx["mean_s"]
        opal["mean_s_rel"] = opal["mean_s"]

    if "ref_z" in opalx:
        opalx["ref_z_rel"] = opalx["ref_z"]
        opal["ref_z_rel"] = opal["ref_z"] - z_offset

    if "t" in opalx:
        opalx["t_rel"] = opalx["t"] - opalx["t"].iloc[0]
        opal["t_rel"] = opal["t"] - opal["t"].iloc[0]

    return opalx, opal


def comparison_table(
    opalx: pd.DataFrame,
    opal: pd.DataFrame,
    *,
    x_column: str,
    columns: list[str],
) -> pd.DataFrame:
    opalx = opalx.sort_values(x_column)
    opal = opal.sort_values(x_column)

    x_min = max(opalx[x_column].min(), opal[x_column].min())
    x_max = min(opalx[x_column].max(), opal[x_column].max())
    mask = (opalx[x_column] >= x_min) & (opalx[x_column] <= x_max)
    base = opalx.loc[mask].copy()

    rows = []
    for column in columns:
        if column not in base or column not in opal:
            continue
        lhs = base[column].to_numpy(dtype=float)
        rhs = np.interp(
            base[x_column].to_numpy(dtype=float),
            opal[x_column].to_numpy(dtype=float),
            opal[column].to_numpy(dtype=float),
        )
        diff = lhs - rhs
        finite = np.isfinite(diff)
        if not finite.any():
            continue
        rows.append(
            {
                "quantity": column,
                "samples": int(finite.sum()),
                "max_abs_diff": float(np.max(np.abs(diff[finite]))),
                "rms_diff": float(np.sqrt(np.mean(diff[finite] ** 2))),
                "final_opalx": float(lhs[finite][-1]),
                "final_opal": float(rhs[finite][-1]),
                "final_diff": float(diff[finite][-1]),
            }
        )

    return pd.DataFrame(rows)


def field_support(df: pd.DataFrame, *, x_column: str, field_column: str, threshold: float) -> tuple[float, float, float]:
    values = df[[x_column, field_column]].to_numpy(dtype=float)
    finite = np.isfinite(values).all(axis=1)
    values = values[finite]
    if values.size == 0:
        return float("nan"), float("nan"), 0.0

    active = np.abs(values[:, 1]) > threshold
    if not active.any():
        return float("nan"), float("nan"), 0.0

    x_active = values[active, 0]
    return float(x_active[0]), float(x_active[-1]), float(x_active[-1] - x_active[0])


def integrate(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def bfield_table(
    opalx: pd.DataFrame,
    opal: pd.DataFrame,
    *,
    x_column: str,
    components: list[tuple[str, str, str]],
    support_fraction: float = 1.0e-6,
) -> pd.DataFrame:
    opalx = opalx.sort_values(x_column)
    opal = opal.sort_values(x_column)

    x_min = max(opalx[x_column].min(), opal[x_column].min())
    x_max = min(opalx[x_column].max(), opal[x_column].max())
    mask = (opalx[x_column] >= x_min) & (opalx[x_column] <= x_max)
    base = opalx.loc[mask].copy()

    rows = []
    for component, opalx_column, opal_column in components:
        if opalx_column not in base or opal_column not in opal:
            continue

        x = base[x_column].to_numpy(dtype=float)
        bx = base[opalx_column].to_numpy(dtype=float)
        bo = np.interp(
            x,
            opal[x_column].to_numpy(dtype=float),
            opal[opal_column].to_numpy(dtype=float),
        )
        finite = np.isfinite(x) & np.isfinite(bx) & np.isfinite(bo)
        if not finite.any():
            continue

        x_finite = x[finite]
        bx = bx[finite]
        bo = bo[finite]
        diff = bx - bo
        shared_peak = max(float(np.max(np.abs(bx))), float(np.max(np.abs(bo))), 1.0)
        threshold = support_fraction * shared_peak
        ox0, ox1, ox_len = field_support(
            opalx, x_column=x_column, field_column=opalx_column, threshold=threshold
        )
        oo0, oo1, oo_len = field_support(
            opal, x_column=x_column, field_column=opal_column, threshold=threshold
        )

        rows.append(
            {
                "component": component,
                "samples": int(finite.sum()),
                "opalx_peak": float(np.max(np.abs(bx))),
                "opal_peak": float(np.max(np.abs(bo))),
                "peak_diff": float(np.max(np.abs(bx)) - np.max(np.abs(bo))),
                "opalx_int": integrate(bx, x_finite),
                "opal_int": integrate(bo, x_finite),
                "int_diff": integrate(diff, x_finite),
                "max_abs_diff": float(np.max(np.abs(diff))),
                "rms_diff": float(np.sqrt(np.mean(diff**2))),
                "opalx_support": ox_len,
                "opal_support": oo_len,
                "opalx_start": ox0,
                "opal_start": oo0,
                "opalx_end": ox1,
                "opal_end": oo1,
            }
        )

    return pd.DataFrame(rows)


def write_plot(
    opalx_design: pd.DataFrame,
    opal_design: pd.DataFrame,
    opalx_stat: pd.DataFrame,
    opal_stat: pd.DataFrame,
    out_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 220,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "font.size": 10,
        }
    )

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(opalx_design["Rz"], opalx_design["Rx"], label="OPALX")
    ax.plot(opal_design["Rz_rel"], opal_design["Rx"], "--", label="OPAL 2022.1")
    ax.set_xlabel("z relative to entrance [m]")
    ax.set_ylabel("x [m]")
    ax.set_title("test-sbend-1 design path")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "design_path_xz.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(opalx_design["s_rel"], opalx_design["Px"], label="OPALX Px")
    ax.plot(opalx_design["s_rel"], opalx_design["Pz"], label="OPALX Pz")
    ax.plot(opal_design["s_rel"], opal_design["Px"], "--", label="OPAL Px")
    ax.plot(opal_design["s_rel"], opal_design["Pz"], "--", label="OPAL Pz")
    ax.set_xlabel("s relative to entrance [m]")
    ax.set_ylabel("normalized momentum")
    ax.set_title("test-sbend-1 design momentum")
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(out_dir / "design_path_momentum.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(opalx_design["s_rel"], opalx_design["Bfy"], label="OPALX Bfy")
    ax.plot(opal_design["s_rel"], opal_design["Bfy"], "--", label="OPAL Bfy")
    ax.set_xlabel("s relative to entrance [m]")
    ax.set_ylabel("sampled vertical B field")
    ax.set_title("test-sbend-1 design-path sampled field")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "design_path_by.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(opalx_stat["s_rel"], opalx_stat["By_ref"], label="OPALX By_ref")
    ax.plot(opal_stat["s_rel"], opal_stat["By_ref"], "--", label="OPAL By_ref")
    ax.set_xlabel("s relative to entrance [m]")
    ax.set_ylabel("sampled vertical B field")
    ax.set_title("test-sbend-1 reference-particle sampled field")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "stat_by_ref.png")
    plt.close(fig)


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No comparable finite quantities found._"

    headers = list(df.columns)
    rows = []
    for _, row in df.iterrows():
        formatted = []
        for header in headers:
            value = row[header]
            if isinstance(value, (float, np.floating)):
                formatted.append(f"{value:.6e}")
            else:
                formatted.append(str(value))
        rows.append(formatted)

    widths = [
        max(len(header), *(len(row[index]) for row in rows))
        for index, header in enumerate(headers)
    ]
    header_line = "| " + " | ".join(
        header.ljust(widths[index]) for index, header in enumerate(headers)
    ) + " |"
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [
        "| " + " | ".join(value.ljust(widths[index]) for index, value in enumerate(row)) + " |"
        for row in rows
    ]
    return "\n".join([header_line, separator, *body])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--opalx-dir", type=Path, required=True)
    parser.add_argument("--opal-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--opal-s-offset", type=float, default=1.0)
    parser.add_argument("--opal-z-offset", type=float, default=1.0)
    parser.add_argument("--plots", action="store_true", help="Also write PNG diagnostic plots.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    opalx_design = read_design_path(args.opalx_dir / "data/test-sbend-1_DesignPath.dat")
    opal_design = read_design_path(args.opal_dir / "data/test-sbend-1_opal_DesignPath.dat")
    opalx_stat = read_sdds_stat(args.opalx_dir / "test-sbend-1.stat")
    opal_stat = read_sdds_stat(args.opal_dir / "test-sbend-1_opal.stat")

    opalx_design, opal_design = add_normalized_coordinates(
        opalx_design, opal_design, s_offset=args.opal_s_offset, z_offset=args.opal_z_offset
    )
    opalx_stat, opal_stat = add_normalized_coordinates(
        opalx_stat, opal_stat, s_offset=args.opal_s_offset, z_offset=args.opal_z_offset
    )

    design_summary = comparison_table(
        opalx_design,
        opal_design,
        x_column="s_rel",
        columns=["Rx", "Ry", "Rz_rel", "Px", "Py", "Pz", "Bfy", "Ekin", "t_rel"],
    )
    stat_summary = comparison_table(
        opalx_stat,
        opal_stat,
        x_column="s_rel",
        columns=[
            "energy",
            "mean_x",
            "mean_y",
            "mean_s_rel",
            "ref_x",
            "ref_y",
            "ref_z_rel",
            "ref_px",
            "ref_py",
            "ref_pz",
            "By_ref",
            "t_rel",
        ],
    )
    design_bfield = bfield_table(
        opalx_design,
        opal_design,
        x_column="s_rel",
        components=[("Bx", "Bfx", "Bfx"), ("By", "Bfy", "Bfy"), ("Bz", "Bfz", "Bfz")],
    )
    stat_bfield = bfield_table(
        opalx_stat,
        opal_stat,
        x_column="s_rel",
        components=[
            ("Bx_ref", "Bx_ref", "Bx_ref"),
            ("By_ref", "By_ref", "By_ref"),
            ("Bz_ref", "Bz_ref", "Bz_ref"),
        ],
    )

    design_summary.to_csv(args.out_dir / "design_path_summary.csv", index=False)
    stat_summary.to_csv(args.out_dir / "stat_summary.csv", index=False)
    design_bfield.to_csv(args.out_dir / "design_path_bfield_summary.csv", index=False)
    stat_bfield.to_csv(args.out_dir / "stat_bfield_summary.csv", index=False)
    if args.plots:
        write_plot(opalx_design, opal_design, opalx_stat, opal_stat, args.out_dir)

    summary = "\n".join(
        [
            "# test-sbend-1 OPALX vs OPAL 2022.1",
            "",
            "Only `*_DesignPath.dat` and `.stat` are compared. HDF5 and load-balance output are ignored.",
            f"OPAL normalization: `s_rel = s - {args.opal_s_offset:g}`, `z_rel = z - {args.opal_z_offset:g}`.",
            "",
            "## Design path",
            "",
            markdown_table(design_summary),
            "",
            "### Sampled design-path B field",
            "",
            "Field peaks are in the output field units and integrals are over `s_rel`.",
            "",
            markdown_table(design_bfield),
            "",
            "## Stat file",
            "",
            markdown_table(stat_summary),
            "",
            "### Sampled reference-particle B field",
            "",
            "This uses the `.stat` `Bx_ref`, `By_ref`, and `Bz_ref` columns.",
            "",
            markdown_table(stat_bfield),
            "",
        ]
    )
    (args.out_dir / "summary.md").write_text(summary)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
