#!/usr/bin/env python3
"""Compare OPALX and old OPAL stat output for test-chicane-distribution-1."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


BEAM_COLUMNS = [
    "energy",
    "dE",
    "rms_x",
    "rms_y",
    "rms_s",
    "rms_px",
    "rms_py",
    "rms_ps",
    "emit_x",
    "emit_y",
    "emit_s",
]

DIAGNOSTIC_COLUMNS = [
    "mean_x",
    "mean_y",
    "mean_s",
    "ref_x",
    "ref_y",
    "ref_z",
    "ref_px",
    "ref_py",
    "ref_pz",
    "By_ref",
]


def read_sdds_stat(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
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
        raise ValueError(f"{path}: could not find SDDS data section")
    if not names:
        raise ValueError(f"{path}: could not find SDDS column definitions")

    rows = []
    for line in lines[data_start:]:
        fields = line.split()
        if len(fields) < len(names):
            continue
        try:
            rows.append([float(value) for value in fields[: len(names)]])
        except ValueError:
            continue

    if not rows:
        raise ValueError(f"{path}: no numeric stat rows found")

    data = np.asarray(rows, dtype=float)
    return {name: data[:, index] for index, name in enumerate(names)}


def comparison_rows(
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray],
        columns: list[str],
        x_column: str = "s") -> list[dict[str, float | int | str]]:
    sx = np.asarray(opalx[x_column], dtype=float)
    so = np.asarray(opal[x_column], dtype=float)
    order_x = np.argsort(sx)
    order_o = np.argsort(so)
    sx = sx[order_x]
    so = so[order_o]

    lo = max(float(np.nanmin(sx)), float(np.nanmin(so)))
    hi = min(float(np.nanmax(sx)), float(np.nanmax(so)))
    mask = (sx >= lo) & (sx <= hi)
    sx_common = sx[mask]

    rows: list[dict[str, float | int | str]] = []
    for column in columns:
        if column not in opalx or column not in opal:
            continue
        vx = np.asarray(opalx[column], dtype=float)[order_x][mask]
        vo = np.asarray(opal[column], dtype=float)[order_o]
        vo_interp = np.interp(sx_common, so, vo)
        diff = vx - vo_interp
        finite = np.isfinite(diff) & np.isfinite(vx) & np.isfinite(vo_interp)
        if not finite.any():
            continue
        rows.append(
            {
                "quantity": column,
                "samples": int(np.count_nonzero(finite)),
                "max_abs_diff": float(np.max(np.abs(diff[finite]))),
                "rms_diff": float(np.sqrt(np.mean(diff[finite] ** 2))),
                "final_opalx": float(vx[finite][-1]),
                "final_opal": float(vo_interp[finite][-1]),
                "final_diff": float(diff[finite][-1]),
            }
        )
    return rows


def field_interval_rows(
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray],
        column: str = "By_ref",
        x_column: str = "s",
        threshold: float = 1.0e-6) -> list[dict[str, float | int | str]]:
    def intervals(data: dict[str, np.ndarray]) -> list[tuple[float, float, int]]:
        s = np.asarray(data[x_column], dtype=float)
        values = np.asarray(data[column], dtype=float)
        mask = np.isfinite(s) & np.isfinite(values) & (np.abs(values) > threshold)
        if not mask.any():
            return []
        starts = np.where(mask & np.r_[True, ~mask[:-1]])[0]
        ends = np.where(mask & np.r_[~mask[1:], True])[0]
        return [(float(s[start]), float(s[end]), int(end - start + 1))
                for start, end in zip(starts, ends)]

    opalx_intervals = intervals(opalx)
    opal_intervals = intervals(opal)
    count = min(len(opalx_intervals), len(opal_intervals))
    rows: list[dict[str, float | int | str]] = []
    for index in range(count):
        x_start, x_end, x_samples = opalx_intervals[index]
        o_start, o_end, o_samples = opal_intervals[index]
        rows.append(
            {
                "interval": index + 1,
                "opalx_start_m": x_start,
                "opal_start_m": o_start,
                "start_diff_mm": 1.0e3 * (x_start - o_start),
                "opalx_end_m": x_end,
                "opal_end_m": o_end,
                "end_diff_mm": 1.0e3 * (x_end - o_end),
                "opalx_samples": x_samples,
                "opal_samples": o_samples,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, float | int | str]]) -> str:
    if not rows:
        return "_No comparable finite quantities found._"
    headers = list(rows[0].keys())
    rendered = []
    for row in rows:
        rendered.append(
            [
                f"{value:.6e}" if isinstance(value, float) else str(value)
                for value in (row[header] for header in headers)
            ]
        )
    widths = [
        max(len(header), *(len(row[index]) for row in rendered))
        for index, header in enumerate(headers)
    ]
    header_line = "| " + " | ".join(
        header.ljust(widths[index]) for index, header in enumerate(headers)
    ) + " |"
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [
        "| " + " | ".join(value.ljust(widths[index]) for index, value in enumerate(row)) + " |"
        for row in rendered
    ]
    return "\n".join([header_line, separator, *body])


def r56_annotation_lines(specs: list[str]) -> list[str]:
    lines = []
    for spec in specs:
        if "=" in spec:
            label, path_text = spec.split("=", 1)
        else:
            path_text = spec
            label = Path(path_text).parent.name
        path = Path(path_text)
        payload = json.loads(path.read_text(encoding="utf-8"))
        lines.append(
            (
                f"{label}: "
                rf"$R_{{56}}={payload['r56_transport_m']:.8e}$ m, "
                rf"$\Delta R_{{56}}={payload['r56_difference_m']:.2e}$ m"
            )
        )
    return lines


def relative_percent(
        opalx_values: np.ndarray,
        opal_values: np.ndarray,
        denominator_floor: float | None = None) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        max_denom = float(np.nanmax(np.abs(opal_values)))
        if denominator_floor is None:
            floor = max(max_denom * 1.0e-14, np.finfo(float).tiny)
        else:
            floor = max(denominator_floor, np.finfo(float).tiny)
        relative = 100.0 * np.abs(opalx_values - opal_values) / np.abs(opal_values)
        relative[np.abs(opal_values) <= floor] = np.nan
    return relative


def positive_limits(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    high = float(np.max(finite))
    if high <= 0.0:
        high = 1.0
    return 0.0, 1.08 * high


def symmetric_limits(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return -1.0, 1.0
    magnitude = max(abs(float(np.min(finite))), abs(float(np.max(finite))))
    if magnitude == 0.0:
        magnitude = 1.0
    return -1.08 * magnitude, 1.08 * magnitude


def draw_line_plot(
        path: Path,
        series: list[tuple[np.ndarray, np.ndarray, str, tuple[int, int, int], bool]],
        xlabel: str,
        ylabel: str,
        title: str,
        y_scale: float = 1.0,
        min_y_span: float | None = None) -> None:
    width, height = 1000, 650
    left, right, top, bottom = 115, 45, 80, 85
    plot_w = width - left - right
    plot_h = height - top - bottom
    font = ImageFont.load_default()
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)

    x_all = np.concatenate([np.asarray(x, dtype=float) for x, _, _, _, _ in series])
    y_all = np.concatenate([np.asarray(y, dtype=float) * y_scale for _, y, _, _, _ in series])
    finite = np.isfinite(x_all) & np.isfinite(y_all)
    x_min, x_max = float(np.min(x_all[finite])), float(np.max(x_all[finite]))
    y_min, y_max = float(np.min(y_all[finite])), float(np.max(y_all[finite]))
    y_span = y_max - y_min
    if min_y_span is not None and y_span < min_y_span:
        center = 0.5 * (y_min + y_max)
        y_min = center - 0.5 * min_y_span
        y_max = center + 0.5 * min_y_span
        y_span = min_y_span
    if y_span == 0.0:
        y_min -= 1.0
        y_max += 1.0
        y_span = 2.0
    y_pad = 0.08 * y_span
    y_min -= y_pad
    y_max += y_pad
    y_span = y_max - y_min

    def sx(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_w

    def sy(value: float) -> float:
        return top + plot_h - (value - y_min) / (y_max - y_min) * plot_h

    def fmt(value: float, span: float | None = None) -> str:
        if span is not None and span < 0.01 * max(1.0, abs(value)):
            return f"{value:.6f}"
        if value == 0.0:
            return "0"
        if abs(value) < 1.0e-2 or abs(value) >= 1.0e3:
            return f"{value:.2e}"
        return f"{value:.4g}"

    for i in range(6):
        xv = x_min + (x_max - x_min) * i / 5.0
        px = sx(xv)
        draw.line((px, top, px, top + plot_h), fill=(225, 225, 225))
        draw.text((px - 20, top + plot_h + 12), fmt(xv), fill=(30, 30, 30), font=font)
        yv = y_min + (y_max - y_min) * i / 5.0
        py = sy(yv)
        draw.line((left, py, left + plot_w, py), fill=(225, 225, 225))
        draw.text((16, py - 6), fmt(yv, y_span), fill=(30, 30, 30), font=font)

    draw.rectangle((left, top, left + plot_w, top + plot_h), outline=(30, 30, 30))

    for x, y, _label, color, dashed in series:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float) * y_scale
        points = [(sx(float(a)), sy(float(b))) for a, b in zip(x, y) if math.isfinite(a + b)]
        if len(points) < 2:
            continue
        if dashed:
            for idx in range(0, len(points) - 1, 3):
                draw.line(points[idx:idx + 2], fill=color, width=2)
        else:
            draw.line(points, fill=color, width=2)

    draw.text((left, 28), title, fill=(0, 0, 0), font=font)
    draw.text((left + plot_w // 2 - 35, height - 35), xlabel, fill=(0, 0, 0), font=font)
    label_img = Image.new("RGBA", (120, 20), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label_img)
    label_draw.text((0, 0), ylabel, fill=(0, 0, 0), font=font)
    label_img = label_img.rotate(90, expand=True)
    img.paste(label_img, (22, top + plot_h // 2 - label_img.height // 2), label_img)

    legend_x = left + plot_w - 195
    legend_y = top + 15
    for idx, (_x, _y, label, color, dashed) in enumerate(series):
        yy = legend_y + idx * 18
        if dashed:
            for xx in range(legend_x, legend_x + 28, 8):
                draw.line((xx, yy + 5, xx + 5, yy + 5), fill=color, width=3)
        else:
            draw.line((legend_x, yy + 5, legend_x + 28, yy + 5), fill=color, width=3)
        draw.text((legend_x + 35, yy), label, fill=(0, 0, 0), font=font)

    img.save(path)


def draw_by_difference_plot(
        path: Path,
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray]) -> None:
    mpl_config_dir = Path(os.environ.get("MPLCONFIGDIR", "/private/tmp/opalx-matplotlib"))
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sx = np.asarray(opalx["s"], dtype=float)
    so = np.asarray(opal["s"], dtype=float)
    bx = np.asarray(opalx["By_ref"], dtype=float)
    bo = np.asarray(opal["By_ref"], dtype=float)

    order_x = np.argsort(sx)
    order_o = np.argsort(so)
    sx, bx = sx[order_x], bx[order_x]
    so, bo = so[order_o], bo[order_o]

    lo = max(float(np.nanmin(sx)), float(np.nanmin(so)))
    hi = min(float(np.nanmax(sx)), float(np.nanmax(so)))
    mask = (sx >= lo) & (sx <= hi)
    sx_common = sx[mask]
    bx_common = bx[mask]
    bo_interp = np.interp(sx_common, so, bo)
    field_floor = max(1.0e-6, 0.5 * float(np.nanmax(np.abs(bo_interp))))
    relative = relative_percent(bx_common, bo_interp, denominator_floor=field_floor)

    with plt.rc_context(
            {
                "font.family": "serif",
                "mathtext.fontset": "stix",
                "font.size": 10,
                "axes.labelsize": 10,
                "axes.titlesize": 11,
                "axes.linewidth": 0.8,
                "xtick.labelsize": 9,
                "ytick.labelsize": 9,
                "legend.fontsize": 8.5,
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "savefig.dpi": 300,
            }):
        fig, ax_field = plt.subplots(figsize=(7.2, 4.2))
        ax_diff = ax_field.twinx()
        fig.subplots_adjust(left=0.12, right=0.86, bottom=0.14, top=0.83)

        field_color = "#009E73"
        diff_color = "#D55E00"
        ax_field.plot(
            sx_common,
            bx_common,
            color=field_color,
            lw=1.65,
            solid_capstyle="round",
        )
        ax_field.plot(
            sx_common,
            bo_interp,
            color=field_color,
            lw=1.25,
            ls=(0, (4.0, 2.0)),
            alpha=0.82,
            dash_capstyle="round",
        )
        ax_diff.plot(
            sx_common,
            relative,
            color=diff_color,
            lw=1.05,
            ls=(0, (1.0, 2.0)),
            alpha=0.9,
            dash_capstyle="round",
        )

        ax_field.set_xlim(lo, hi)
        ax_diff.set_ylim(*positive_limits(relative))
        ax_field.set_xlabel(r"$s$ [m]")
        ax_field.set_ylabel(r"$B_{y,\mathrm{ref}}$ [T]")
        ax_diff.set_ylabel(r"$|\Delta|$ [%]")
        ax_diff.yaxis.set_label_coords(1.10, 0.5)
        ax_diff.yaxis.label.set_color("#4c4c4c")
        ax_diff.tick_params(axis="y", colors="#4c4c4c")

        for axis in (ax_field, ax_diff):
            axis.grid(True, which="major", color="#d9d9d9", lw=0.6)
            axis.grid(True, which="minor", color="#eeeeee", lw=0.4)
            axis.minorticks_on()
            axis.spines["top"].set_visible(False)
        ax_diff.axhline(0.0, color="#5f5f5f", lw=0.6, alpha=0.55, zorder=1)
        ax_field.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)
        ax_diff.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)

        handles = [
            Line2D([0], [0], color=field_color, lw=1.65, label=r"OPALX $B_{y,\mathrm{ref}}$"),
            Line2D(
                [0],
                [0],
                color=field_color,
                lw=1.25,
                ls=(0, (4.0, 2.0)),
                label=r"OPAL $B_{y,\mathrm{ref}}$",
            ),
            Line2D([0], [0], color=diff_color, lw=1.05, ls=(0, (1.0, 2.0)), label=r"$|\Delta|$ [%]"),
        ]
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.96),
            ncol=3,
            frameon=False,
            borderpad=0.0,
            columnspacing=1.7,
            handlelength=2.4,
        )

        fig.savefig(path, bbox_inches="tight")
        if path.suffix.lower() == ".png":
            fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def draw_reference_orbit_difference_plot(
        path: Path,
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray]) -> None:
    mpl_config_dir = Path(os.environ.get("MPLCONFIGDIR", "/private/tmp/opalx-matplotlib"))
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sx = np.asarray(opalx["s"], dtype=float)
    so = np.asarray(opal["s"], dtype=float)
    zx = np.asarray(opalx["ref_z"], dtype=float)
    zo = np.asarray(opal["ref_z"], dtype=float)
    xx = np.asarray(opalx["ref_x"], dtype=float)
    xo = np.asarray(opal["ref_x"], dtype=float)

    order_x = np.argsort(sx)
    order_o = np.argsort(so)
    sx, zx, xx = sx[order_x], zx[order_x], xx[order_x]
    so, zo, xo = so[order_o], zo[order_o], xo[order_o]

    lo = max(float(np.nanmin(sx)), float(np.nanmin(so)))
    hi = min(float(np.nanmax(sx)), float(np.nanmax(so)))
    mask = (sx >= lo) & (sx <= hi)
    sx_common = sx[mask]
    zx_common = zx[mask]
    xx_common = xx[mask]
    zo_interp = np.interp(sx_common, so, zo)
    xo_interp = np.interp(sx_common, so, xo)
    dx_mm = 1.0e3 * (xx_common - xo_interp)
    dz_mm = 1.0e3 * (zx_common - zo_interp)

    with plt.rc_context(
            {
                "font.family": "serif",
                "mathtext.fontset": "stix",
                "font.size": 10,
                "axes.labelsize": 10,
                "axes.titlesize": 11,
                "axes.linewidth": 0.8,
                "xtick.labelsize": 9,
                "ytick.labelsize": 9,
                "legend.fontsize": 8.5,
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "savefig.dpi": 300,
            }):
        fig, (ax_orbit, ax_diff) = plt.subplots(
            2,
            1,
            figsize=(7.2, 5.25),
            gridspec_kw={"height_ratios": [2.2, 1.0]},
        )
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.11, top=0.86, hspace=0.32)

        orbit_color = "#0072B2"
        dx_color = "#0072B2"
        dz_color = "#D55E00"
        ax_orbit.plot(
            zx,
            1.0e3 * xx,
            color=orbit_color,
            lw=1.8,
            solid_capstyle="round",
        )
        ax_orbit.plot(
            zo,
            1.0e3 * xo,
            color=orbit_color,
            lw=1.35,
            ls=(0, (4.0, 2.0)),
            alpha=0.82,
            dash_capstyle="round",
        )

        ax_diff.plot(sx_common, dx_mm, color=dx_color, lw=1.25, label=r"$\Delta x_{\mathrm{ref}}$")
        ax_diff.plot(
            sx_common,
            dz_mm,
            color=dz_color,
            lw=1.15,
            ls=(0, (4.0, 2.0)),
            label=r"$\Delta z_{\mathrm{ref}}$",
        )
        ax_diff.axhline(0.0, color="#5f5f5f", lw=0.65, alpha=0.65, zorder=1)

        ax_orbit.set_xlabel(r"$z_{\mathrm{ref}}$ [m]")
        ax_orbit.set_ylabel(r"$x_{\mathrm{ref}}$ [mm]")
        ax_diff.set_xlabel(r"$s$ [m]")
        ax_diff.set_ylabel(r"OPALX $-$ OPAL [mm]")

        for axis in (ax_orbit, ax_diff):
            axis.grid(True, which="major", color="#d9d9d9", lw=0.6)
            axis.grid(True, which="minor", color="#eeeeee", lw=0.4)
            axis.minorticks_on()
            axis.spines["top"].set_visible(False)
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)

        ax_diff.set_xlim(lo, hi)
        ax_diff.set_ylim(*symmetric_limits(np.concatenate((dx_mm, dz_mm))))

        handles = [
            Line2D([0], [0], color=orbit_color, lw=1.8, label="OPALX"),
            Line2D([0], [0], color=orbit_color, lw=1.35, ls=(0, (4.0, 2.0)), label="OPAL"),
            Line2D([0], [0], color=dx_color, lw=1.25, label=r"$\Delta x_{\mathrm{ref}}$"),
            Line2D(
                [0],
                [0],
                color=dz_color,
                lw=1.15,
                ls=(0, (4.0, 2.0)),
                label=r"$\Delta z_{\mathrm{ref}}$",
            ),
        ]
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.98),
            ncol=4,
            frameon=False,
            columnspacing=1.35,
            handlelength=2.2,
        )

        fig.savefig(path, bbox_inches="tight")
        if path.suffix.lower() == ".png":
            fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def draw_rms_relative_difference_plot(
        path: Path,
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray],
        r56_lines: list[str] | None = None) -> None:
    mpl_config_dir = Path(os.environ.get("MPLCONFIGDIR", "/private/tmp/opalx-matplotlib"))
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sx = np.asarray(opalx["s"], dtype=float)
    so = np.asarray(opal["s"], dtype=float)
    order_x = np.argsort(sx)
    order_o = np.argsort(so)
    sx = sx[order_x]
    so = so[order_o]

    lo = max(float(np.nanmin(sx)), float(np.nanmin(so)))
    hi = min(float(np.nanmax(sx)), float(np.nanmax(so)))
    mask = (sx >= lo) & (sx <= hi)
    sx_common = sx[mask]

    quantities = [
        ("rms_x", r"$\sigma_x$", "#0072B2"),
        ("rms_y", r"$\sigma_y$", "#009E73"),
        ("rms_s", r"$\sigma_s$", "#D55E00"),
    ]
    series: list[tuple[str, str, str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for quantity, label, color in quantities:
        if quantity not in opalx or quantity not in opal:
            continue
        opalx_values = np.asarray(opalx[quantity], dtype=float)[order_x][mask]
        opal_values = np.asarray(opal[quantity], dtype=float)[order_o]
        opal_interp = np.interp(sx_common, so, opal_values)
        opalx_mm = 1.0e3 * opalx_values
        opal_mm = 1.0e3 * opal_interp
        absolute_mm = opalx_mm - opal_mm
        relative = relative_percent(opalx_values, opal_interp)
        series.append((quantity, label, color, opalx_mm, opal_mm, absolute_mm, relative))

    if not series:
        raise ValueError("No RMS columns available for OPALX/OPAL comparison")

    with plt.rc_context(
            {
                "font.family": "serif",
                "mathtext.fontset": "stix",
                "font.size": 10,
                "axes.labelsize": 10,
                "axes.titlesize": 11,
                "axes.linewidth": 0.8,
                "xtick.labelsize": 9,
                "ytick.labelsize": 9,
                "legend.fontsize": 8.5,
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "savefig.dpi": 300,
            }):
        fig, (ax_values, ax_rel) = plt.subplots(
            2,
            1,
            figsize=(7.2, 5.25),
            sharex=True,
            gridspec_kw={"height_ratios": [2.4, 1.0]},
        )
        fig.subplots_adjust(left=0.11, right=0.88, bottom=0.12, top=0.73 if r56_lines else 0.81, hspace=0.08)

        for _quantity, label, color, opalx_values, opal_values, absolute, relative in series:
            ax_values.plot(
                sx_common,
                opalx_values,
                color=color,
                lw=1.8,
                solid_capstyle="round",
                label=f"OPALX {label}",
            )
            ax_values.plot(
                sx_common,
                opal_values,
                color=color,
                lw=1.35,
                ls=(0, (4.0, 2.0)),
                alpha=0.82,
                dash_capstyle="round",
                label=f"OPAL {label}",
            )
            ax_rel.plot(
                sx_common,
                relative,
                color=color,
                lw=1.1,
                ls=(0, (1.0, 2.0)),
                alpha=0.9,
                dash_capstyle="round",
                label=f"{label} rel.",
            )

        rel_values = np.concatenate([relative for *_unused, relative in series])
        ax_rel.set_ylim(*positive_limits(rel_values))

        ax_values.set_xlim(lo, hi)
        ax_rel.set_xlabel(r"$s$ [m]")
        ax_values.set_ylabel("RMS beam size [mm]")
        ax_rel.set_ylabel(r"$|\Delta|$ [%]")

        for axis in (ax_values, ax_rel):
            axis.grid(True, which="major", color="#d9d9d9", lw=0.6)
            axis.grid(True, which="minor", color="#eeeeee", lw=0.4)
            axis.minorticks_on()
        ax_rel.axhline(0.0, color="#5f5f5f", lw=0.65, alpha=0.6, zorder=1)
        ax_values.spines["top"].set_visible(False)
        ax_rel.spines["top"].set_visible(False)
        ax_rel.tick_params(axis="y", colors="#4c4c4c")
        ax_rel.yaxis.tick_right()
        ax_rel.yaxis.set_label_position("right")
        ax_rel.yaxis.set_label_coords(1.09, 0.5)
        ax_rel.yaxis.label.set_color("#4c4c4c")
        ax_values.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)
        ax_rel.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)

        quantity_handles = [
            Line2D([0], [0], color=color, lw=2.0, label=label)
            for _quantity, label, color in quantities
        ]
        style_handles = [
            Line2D([0], [0], color="#333333", lw=1.8, label="OPALX"),
            Line2D([0], [0], color="#333333", lw=1.35, ls=(0, (4.0, 2.0)), label="OPAL"),
            Line2D([0], [0], color="#333333", lw=1.1, ls=(0, (1.0, 2.0)), label=r"$|\Delta|$ [%]"),
        ]
        legend = fig.legend(
            handles=[*quantity_handles, *style_handles],
            loc="upper center",
            bbox_to_anchor=(0.5, 0.98),
            ncol=6,
            frameon=False,
            borderpad=0.0,
            columnspacing=1.45,
            handlelength=2.4,
        )
        for line in legend.get_lines():
            line.set_linewidth(1.8)
        if r56_lines:
            fig.text(
                0.5,
                0.89,
                "\n".join(r56_lines),
                ha="center",
                va="top",
                fontsize=8.5,
            )

        fig.savefig(path, dpi=300, bbox_inches="tight")
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def draw_three_panel_difference_plot(
        path: Path,
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray],
        quantities: list[tuple[str, str, str]],
        left_ylabel: str,
        right_ylabel: str) -> None:
    mpl_config_dir = Path(os.environ.get("MPLCONFIGDIR", "/private/tmp/opalx-matplotlib"))
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sx = np.asarray(opalx["s"], dtype=float)
    so = np.asarray(opal["s"], dtype=float)
    order_x = np.argsort(sx)
    order_o = np.argsort(so)
    sx = sx[order_x]
    so = so[order_o]

    lo = max(float(np.nanmin(sx)), float(np.nanmin(so)))
    hi = min(float(np.nanmax(sx)), float(np.nanmax(so)))
    mask = (sx >= lo) & (sx <= hi)
    sx_common = sx[mask]

    series = []
    for quantity, label, color in quantities:
        if quantity not in opalx or quantity not in opal:
            continue
        opalx_values = np.asarray(opalx[quantity], dtype=float)[order_x][mask]
        opal_values = np.asarray(opal[quantity], dtype=float)[order_o]
        opal_interp = np.interp(sx_common, so, opal_values)
        difference = opalx_values - opal_interp
        relative = relative_percent(opalx_values, opal_interp)
        series.append((quantity, label, color, opalx_values, opal_interp, difference, relative))

    if not series:
        raise ValueError("No requested stat columns available for OPALX/OPAL comparison")

    with plt.rc_context(
            {
                "font.family": "serif",
                "mathtext.fontset": "stix",
                "font.size": 9.5,
                "axes.labelsize": 9.5,
                "axes.titlesize": 10,
                "axes.linewidth": 0.8,
                "xtick.labelsize": 8.5,
                "ytick.labelsize": 8.5,
                "legend.fontsize": 8.5,
                "pdf.fonttype": 42,
                "ps.fonttype": 42,
                "savefig.dpi": 300,
            }):
        fig, axes = plt.subplots(
            len(series),
            2,
            figsize=(8.4, 5.8),
            sharex=True,
            gridspec_kw={"width_ratios": [1.35, 1.0]},
        )
        if len(series) == 1:
            axes = np.asarray([axes])
        fig.subplots_adjust(left=0.10, right=0.87, bottom=0.12, top=0.87, hspace=0.12, wspace=0.30)

        for row, (_quantity, label, color, opalx_values, opal_values, difference, relative) in enumerate(series):
            ax_value = axes[row, 0]
            ax_rel = axes[row, 1]
            ax_value.plot(sx_common, opalx_values, color=color, lw=1.65, solid_capstyle="round")
            ax_value.plot(
                sx_common,
                opal_values,
                color=color,
                lw=1.25,
                ls=(0, (4.0, 2.0)),
                alpha=0.82,
                dash_capstyle="round",
            )
            ax_rel.plot(
                sx_common,
                relative,
                color=color,
                lw=1.0,
                ls=(0, (1.0, 2.0)),
                alpha=0.9,
                dash_capstyle="round",
            )

            ax_rel.set_ylim(*positive_limits(relative))

            for axis in (ax_value, ax_rel):
                axis.set_xlim(lo, hi)
                axis.grid(True, which="major", color="#d9d9d9", lw=0.6)
                axis.grid(True, which="minor", color="#eeeeee", lw=0.4)
                axis.minorticks_on()
            ax_rel.axhline(0.0, color="#5f5f5f", lw=0.6, alpha=0.55, zorder=1)
            ax_value.spines["top"].set_visible(False)
            ax_rel.spines["top"].set_visible(False)
            ax_rel.tick_params(axis="y", colors="#4c4c4c")
            ax_rel.yaxis.tick_right()
            ax_rel.yaxis.set_label_position("right")
            ax_value.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)
            ax_rel.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True)
            ax_value.text(
                0.012,
                0.84,
                label,
                transform=ax_value.transAxes,
                color=color,
                fontsize=10,
                ha="left",
                va="top",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5},
            )
            ax_rel.text(
                0.02,
                0.84,
                label,
                transform=ax_rel.transAxes,
                color=color,
                fontsize=10,
                ha="left",
                va="top",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5},
            )

        axes[-1, 0].set_xlabel(r"$s$ [m]")
        axes[-1, 1].set_xlabel(r"$s$ [m]")
        fig.text(0.03, 0.5, left_ylabel, rotation=90, ha="center", va="center")
        fig.text(
            0.965,
            0.5,
            right_ylabel,
            rotation=90,
            ha="center",
            va="center",
            color="#4c4c4c",
        )

        legend_handles = [
            Line2D([0], [0], color="#333333", lw=1.65, label="OPALX"),
            Line2D([0], [0], color="#333333", lw=1.25, ls=(0, (4.0, 2.0)), label="OPAL"),
            Line2D([0], [0], color="#333333", lw=1.0, ls=(0, (1.0, 2.0)), label=r"$|\Delta|$ [%]"),
        ]
        legend = fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=3,
            frameon=False,
            borderpad=0.0,
            columnspacing=1.8,
            handlelength=2.5,
        )
        for line in legend.get_lines():
            line.set_linewidth(1.8)

        fig.savefig(path, dpi=300, bbox_inches="tight")
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def write_plots(
        opalx: dict[str, np.ndarray],
        opal: dict[str, np.ndarray],
        out_dir: Path,
        r56_lines: list[str] | None = None,
        beam_only: bool = False) -> None:
    colors = {
        "x": (31, 119, 180),
        "y": (44, 160, 44),
        "z": (148, 103, 189),
        "s": (214, 39, 40),
        "energy": (127, 127, 127),
    }
    draw_line_plot(
        out_dir / "stat_rms_opalx_opal.png",
        [
            (opalx["s"], opalx["rms_x"], "OPALX rms_x", colors["x"], False),
            (opal["s"], opal["rms_x"], "OPAL rms_x", colors["x"], True),
            (opalx["s"], opalx["rms_y"], "OPALX rms_y", colors["y"], False),
            (opal["s"], opal["rms_y"], "OPAL rms_y", colors["y"], True),
            (opalx["s"], opalx["rms_s"], "OPALX rms_s", colors["s"], False),
            (opal["s"], opal["rms_s"], "OPAL rms_s", colors["s"], True),
        ],
        "s [m]",
        "RMS [mm]",
        "OPALX vs OPAL stat RMS sizes",
        y_scale=1.0e3,
    )
    draw_rms_relative_difference_plot(
        out_dir / "stat_rms_relative_opalx_opal.png",
        opalx,
        opal,
        r56_lines,
    )
    draw_three_panel_difference_plot(
        out_dir / "stat_rms_momentum_difference_opalx_opal.png",
        opalx,
        opal,
        [
            ("rms_px", r"$\sigma_{p_x}$", "#0072B2"),
            ("rms_py", r"$\sigma_{p_y}$", "#009E73"),
            ("rms_ps", r"$\sigma_{p_s}$", "#D55E00"),
        ],
        r"RMS momentum [$\beta\gamma$]",
        r"$|\Delta|$ [%]",
    )
    draw_three_panel_difference_plot(
        out_dir / "stat_emittance_difference_opalx_opal.png",
        opalx,
        opal,
        [
            ("emit_x", r"$\varepsilon_x$", "#0072B2"),
            ("emit_y", r"$\varepsilon_y$", "#009E73"),
            ("emit_s", r"$\varepsilon_s$", "#D55E00"),
        ],
        "RMS emittance",
        r"$|\Delta|$ [%]",
    )
    draw_line_plot(
        out_dir / "stat_energy_opalx_opal.png",
        [
            (opalx["s"], opalx["energy"], "OPALX energy", colors["energy"], False),
            (opal["s"], opal["energy"], "OPAL energy", colors["energy"], True),
        ],
        "s [m]",
        "energy [MeV]",
        "OPALX vs OPAL stat energy",
        min_y_span=1.0e-3,
    )
    if "dE" in opalx and "dE" in opal:
        draw_line_plot(
            out_dir / "stat_dE_opalx_opal.png",
            [
                (opalx["s"], opalx["dE"], "OPALX dE", colors["energy"], False),
                (opal["s"], opal["dE"], "OPAL dE", colors["energy"], True),
            ],
            "s [m]",
            "dE [MeV]",
            "OPALX vs OPAL stat energy spread",
        )
    if beam_only:
        return
    if "ref_x" in opalx and "ref_z" in opalx and "ref_x" in opal and "ref_z" in opal:
        draw_reference_orbit_difference_plot(
            out_dir / "stat_reference_xz_opalx_opal.png",
            opalx,
            opal,
        )
    if "By_ref" in opalx and "By_ref" in opal:
        draw_by_difference_plot(out_dir / "stat_bfield_ref_opalx_opal.png", opalx, opal)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opalx-stat", type=Path, default=script_dir / "test-chicane-distribution-1.stat")
    parser.add_argument("--opal-stat", type=Path, default=script_dir / "test-chicane-distribution-1_opal.stat")
    parser.add_argument("--out-dir", type=Path, default=script_dir / "comparison" / "opal_code_compare")
    parser.add_argument(
        "--beam-only",
        action="store_true",
        help="Compare only beam quantities from the stat file: energy, dE, RMS positions/momenta, and emittances.",
    )
    parser.add_argument(
        "--r56-summary",
        action="append",
        default=[],
        help="Optional LABEL=path to comparison_summary.json written by compare_results.py.",
    )
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    opalx = read_sdds_stat(args.opalx_stat)
    opal = read_sdds_stat(args.opal_stat)

    columns = BEAM_COLUMNS if args.beam_only else [*BEAM_COLUMNS, *DIAGNOSTIC_COLUMNS]
    rows = comparison_rows(opalx, opal, columns)
    by_interval_rows = [] if args.beam_only else field_interval_rows(opalx, opal)
    write_csv(args.out_dir / "stat_summary.csv", rows)
    if not args.beam_only:
        write_csv(args.out_dir / "by_field_intervals.csv", by_interval_rows)
    r56_lines = [] if args.beam_only else r56_annotation_lines(args.r56_summary)
    if args.plots:
        write_plots(opalx, opal, args.out_dir, r56_lines, beam_only=args.beam_only)

    summary_lines = [
        "# test-chicane-distribution-1 OPALX vs OPAL 2022.1",
        "",
        f"OPALX stat: `{args.opalx_stat}`",
        f"OPAL stat: `{args.opal_stat}`",
        "",
        "Both files are compared on the shared `s` range by interpolating OPAL onto OPALX samples.",
        "",
    ]
    if args.beam_only:
        summary_lines.extend(
            [
                "Comparison mode: stat beam quantities only. Temporal monitors, reference-particle fields, and R56 fits are not used.",
                "",
            ]
        )
    summary_lines.extend([markdown_table(rows), ""])
    if not args.beam_only:
        summary_lines.extend(
            [
            "## By Field Intervals",
            "",
            (
                "Intervals where `abs(By_ref) > 1e-6 T`; differences are `OPALX - OPAL`. "
                "A one-step shift is about `2.99792458 mm` for `DT = 1e-11 s`."
            ),
            "",
            markdown_table(by_interval_rows),
            "",
            ]
        )
    plot_lines = [
        "- `stat_rms_opalx_opal.png`",
        "- `stat_rms_relative_opalx_opal.png`",
        "- `stat_rms_relative_opalx_opal.pdf`",
        "- `stat_rms_momentum_difference_opalx_opal.png`",
        "- `stat_rms_momentum_difference_opalx_opal.pdf`",
        "- `stat_emittance_difference_opalx_opal.png`",
        "- `stat_emittance_difference_opalx_opal.pdf`",
        "- `stat_energy_opalx_opal.png`",
    ]
    if "dE" in opalx and "dE" in opal:
        plot_lines.append("- `stat_dE_opalx_opal.png`")
    if not args.beam_only:
        plot_lines.extend(
            [
                "- `stat_reference_xz_opalx_opal.png`",
                "- `stat_bfield_ref_opalx_opal.png`",
                "- `stat_bfield_ref_opalx_opal.pdf`",
            ]
        )
    summary_lines.extend(
        [
            "## Plots",
            "",
            *plot_lines,
            "",
        ]
    )
    summary = "\n".join(summary_lines)
    if r56_lines:
        summary += "\n## R56 Fits\n\n" + "\n".join(f"- {line}" for line in r56_lines) + "\n"
    (args.out_dir / "summary.md").write_text(summary + "\n", encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
