#!/usr/bin/env python3
"""Build a compact PDF overview of the sandbox regression comparisons.

Defaults write the standalone overview TeX/PDF, metrics, baseline, and figures
under ``sandbox/note`` / ``sandbox/note/figs``.  This keeps the overview tooling
usable after the old ``sandbox/regression`` directory is removed.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
NOTE_DIR = ROOT / "sandbox/note"
NOTE_FIG_DIR = NOTE_DIR / "figs"
DEFAULT_CSV = NOTE_DIR / "current_metrics.csv"
DEFAULT_TEX = NOTE_DIR / "sandbox_regression_overview.tex"
DEFAULT_BASELINE = NOTE_DIR / "sandbox_regression_baseline.json"
DEFAULT_IMPACT_PLOT = NOTE_FIG_DIR / "opalx_impact_drift_comparison.png"
DEFAULT_LOSS_PLOT = NOTE_FIG_DIR / "gamma_gamma_cylinder_losses.png"
DEFAULT_LARGE_CYLINDER_PLOT = NOTE_FIG_DIR / "gamma_gamma_large_cylinder_crossings.png"
DEFAULT_LARGE_CYLINDER_CHARGE_COMPARE_PLOT = NOTE_FIG_DIR / "gamma_gamma_large_cylinder_charge_compare.png"


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
        "<": r"\textless{}",
        ">": r"\textgreater{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def metric_map(csv_path: Path) -> dict[str, object]:
    frame = pd.read_csv(csv_path)
    if frame.empty:
        raise ValueError(f"{csv_path} is empty")
    return dict(zip(frame["metric"], frame["value"], strict=True))


def baseline_metadata(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    metadata = data.get("_metadata", {})
    if not isinstance(metadata, dict):
        return {}
    return metadata


def read_sdds_stat(path: Path) -> tuple[pd.DataFrame, list[str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    lines = path.read_text(errors="replace").splitlines()

    columns: list[str] = []
    data_start: int | None = None
    idx = 0
    while idx < len(lines):
        tag = lines[idx].strip()
        if tag == "&column":
            idx += 1
            name: str | None = None
            while idx < len(lines) and lines[idx].strip() != "&end":
                stripped = lines[idx].strip().rstrip(",")
                if stripped.startswith("name="):
                    name = stripped.split("=", 1)[1].strip().strip('"')
                idx += 1
            if name is not None:
                columns.append(name)
        elif tag == "&data":
            while idx < len(lines) and lines[idx].strip() != "&end":
                idx += 1
            data_start = idx + 1
            break
        idx += 1

    if data_start is None:
        raise ValueError(f"missing &data block in {path}")
    if data_start + 3 > len(lines):
        raise ValueError(f"missing SDDS parameter rows in {path}")

    data_lines = [line for line in lines[data_start + 3:] if line.strip()]
    rows = [[float(value) for value in line.split()] for line in data_lines]
    return pd.DataFrame(rows, columns=columns), columns


def make_impact_plot(output_path: Path) -> None:
    env = {
        **os.environ,
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(ROOT / "sandbox/data/.plot-cache/matplotlib"),
        "XDG_CACHE_HOME": str(ROOT / "sandbox/data/.plot-cache/xdg"),
    }
    subprocess.run(
        [
            str(ROOT / ".venv-h6/bin/python"),
            str(ROOT / "sandbox/OPALX-IMPACT/plot_extracted_graph.py"),
            "--input",
            str(ROOT / "sandbox/OPALX-IMPACT/extracted_graph_data.csv"),
            "--input-opalx",
            str(ROOT / "sandbox/OPALX-IMPACT/drift-1.stat"),
            "--output",
            str(output_path),
        ],
        cwd=ROOT,
        env=env,
        check=True,
    )


def make_cylinder_crossing_plot(
        output_path: Path,
        stem: str,
        aperture_radius_m: float,
        cylinder_begin_m: float,
        cylinder_end_m: float,
        x_range_m: tuple[float, float],
        bins: int,
) -> None:
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "sandbox/data/.plot-cache/matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "sandbox/data/.plot-cache/xdg"))

    import h5py
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    import numpy as np

    track_dir = ROOT / "sandbox/track-e-p"

    fig, (ax, cyl) = plt.subplots(
        2,
        1,
        figsize=(6.6, 3.6),
        height_ratios=(2.0, 0.72),
        sharex=True,
        constrained_layout=True,
    )
    species = [
        (f"{stem}_c1.h5", r"$e^-$", "#2AA6B8"),
        (f"{stem}_c2.h5", r"$e^+$", "#D33682"),
    ]
    max_crossings = 0.0
    for filename, label, color in species:
        first_crossing_z: dict[int, float] = {}
        with h5py.File(track_dir / filename, "r") as h5:
            steps = sorted(
                    (key for key in h5.keys() if key.startswith("Step#")),
                    key=lambda key: int(key.split("#", 1)[1]))
            for step in steps:
                group = h5[step]
                ids = np.asarray(group["id"], dtype=int)
                x = np.asarray(group["x"], dtype=float)
                y = np.asarray(group["y"], dtype=float)
                z = np.asarray(group["z"], dtype=float)
                outside = np.hypot(x, y) >= aperture_radius_m
                for particle_id, z_position in zip(ids[outside], z[outside], strict=True):
                    first_crossing_z.setdefault(int(particle_id), float(z_position))

        if first_crossing_z:
            crossing_z = np.fromiter(first_crossing_z.values(), dtype=float)
            edges = np.linspace(x_range_m[0], x_range_m[1], bins + 1)
            counts, edges = np.histogram(crossing_z, bins=edges)
            centers = 0.5 * (edges[:-1] + edges[1:])
            width = 0.42 * (edges[1] - edges[0])
            offset = -0.5 * width if label == r"$e^-$" else 0.5 * width
            max_crossings = max(max_crossings, float(counts.max(initial=0)))
        else:
            centers = np.array([])
            counts = np.array([])
            width = 0.0016
            offset = 0.0
        ax.bar(
            centers + offset,
            counts,
            width=width,
            align="center",
            alpha=0.78,
            color=color,
            edgecolor="black",
            linewidth=0.35,
            label=label,
        )

    ax.set_ylabel("first crossings")
    ax.legend(loc="upper right", frameon=False, ncols=2)
    ax.grid(axis="y", color="0.86", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    if max_crossings == 0.0:
        ax.text(
            0.5,
            0.52,
            "no cylinder-edge crossings recorded",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="0.35",
        )
        ax.set_ylim(0.0, 1.0)

    z_min, z_max = cylinder_begin_m, cylinder_end_m
    cyl.set_ylim(-1.0, 1.0)
    cyl.set_yticks([])
    cyl.set_xlabel("z [m]")
    cyl.spines[["top", "right", "left"]].set_visible(False)
    cyl.plot([z_min, z_max], [0.62, 0.62], color="black", linewidth=1.0)
    cyl.plot([z_min, z_max], [-0.62, -0.62], color="black", linewidth=1.0)
    left = Ellipse(
            (z_min, 0.0), width=0.0045, height=1.24, fill=False, color="black",
            linewidth=1.0)
    right = Ellipse(
            (z_max, 0.0), width=0.0045, height=1.24, fill=False, color="black",
            linewidth=1.0)
    cyl.add_patch(left)
    cyl.add_patch(right)
    ip = 0.5 * (cylinder_begin_m + cylinder_end_m)
    cyl.axvline(ip, color="0.2", linewidth=0.8)
    cyl.text(ip, -0.86, "IP", ha="center", va="top")
    cyl.text(ip, 0.77, "cylindrical BeamBeam aperture", ha="center", va="bottom", fontsize=8)
    ax.set_xlim(*x_range_m)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def first_crossing_positions(stem: str, container: str, aperture_radius_m: float) -> pd.Series:
    import h5py
    import numpy as np

    path = ROOT / "sandbox/track-e-p" / f"{stem}_{container}.h5"
    first_crossing_z: dict[int, float] = {}
    with h5py.File(path, "r") as h5:
        steps = sorted(
                (key for key in h5.keys() if key.startswith("Step#")),
                key=lambda key: int(key.split("#", 1)[1]))
        for step in steps:
            group = h5[step]
            ids = np.asarray(group["id"], dtype=int)
            x = np.asarray(group["x"], dtype=float)
            y = np.asarray(group["y"], dtype=float)
            z = np.asarray(group["z"], dtype=float)
            outside = np.hypot(x, y) >= aperture_radius_m
            for particle_id, z_position in zip(ids[outside], z[outside], strict=True):
                first_crossing_z.setdefault(int(particle_id), float(z_position))
    return pd.Series(first_crossing_z, dtype=float).sort_index()


def make_cylinder_loss_plot(output_path: Path) -> None:
    make_cylinder_crossing_plot(
            output_path,
            stem="gamma_gamma_pairs-2",
            aperture_radius_m=0.01,
            cylinder_begin_m=0.01,
            cylinder_end_m=0.06,
            x_range_m=(0.0, 0.07),
            bins=14)


def make_large_cylinder_crossing_plot(output_path: Path) -> None:
    make_cylinder_crossing_plot(
            output_path,
            stem="gamma_gamma_pairs-large-cylinder",
            aperture_radius_m=0.15,
            cylinder_begin_m=0.01,
            cylinder_end_m=0.33,
            x_range_m=(-0.06, 0.47),
            bins=18)


def make_large_cylinder_charge_compare_plot(output_path: Path) -> None:
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "sandbox/data/.plot-cache/matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "sandbox/data/.plot-cache/xdg"))

    import matplotlib.pyplot as plt
    import numpy as np

    nominal_stem = "gamma_gamma_pairs-large-cylinder"
    reduced_stem = "gamma_gamma_pairs-large-cylinder-q1em5"
    edges = np.linspace(-0.06, 0.47, 19)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]

    fig, axes = plt.subplots(2, 1, figsize=(6.8, 4.5), sharex=True, constrained_layout=True)
    summary_rows: list[tuple[str, int, int, int, float]] = []
    for ax, container, label, color in [
        (axes[0], "c1", r"$e^-$", "#2AA6B8"),
        (axes[1], "c2", r"$e^+$", "#D33682"),
    ]:
        nominal = first_crossing_positions(nominal_stem, container, 0.15)
        reduced = first_crossing_positions(reduced_stem, container, 0.15)
        nominal_counts, _ = np.histogram(nominal.to_numpy(), bins=edges)
        reduced_counts, _ = np.histogram(reduced.to_numpy(), bins=edges)
        diff = reduced_counts - nominal_counts
        shared_ids = nominal.index.intersection(reduced.index)
        max_abs_shift_um = 0.0
        if len(shared_ids) > 0:
            max_abs_shift_um = float((reduced.loc[shared_ids] - nominal.loc[shared_ids]).abs().max() * 1.0e6)
        summary_rows.append(
                (container, int(nominal_counts.sum()), int(reduced_counts.sum()),
                 int(np.abs(diff).sum()), max_abs_shift_um))

        ax.step(centers, nominal_counts, where="mid", color=color, linewidth=1.8, label="nominal charge")
        ax.step(
                centers,
                reduced_counts,
                where="mid",
                color="black",
                linewidth=1.2,
                linestyle="--",
                label=r"$10^{-5}$ primary charge")
        ax.bar(
                centers,
                diff,
                width=0.58 * width,
                color="0.72",
                edgecolor="0.35",
                linewidth=0.25,
                alpha=0.75,
                label="reduced - nominal")
        ax.axhline(0.0, color="0.25", linewidth=0.6)
        ax.set_ylabel(f"{label}\ncounts")
        ax.grid(axis="y", color="0.88", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].legend(loc="upper right", frameon=False, ncols=3, fontsize=8)
    axes[1].set_xlabel("first cylinder-edge crossing z [m]")
    fig.suptitle("15 cm cylinder crossing histogram comparison")
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def fmt(metrics: dict[str, object], key: str, unit: str = "") -> str:
    value = metrics[key]
    if isinstance(value, str):
        return latex_escape(value)
    number = float(value)
    if math.isnan(number):
        rendered = "nan"
    elif number == 0:
        rendered = "0"
    elif 1.0e-3 <= abs(number) < 1.0e4:
        rendered = f"{number:.6g}"
    else:
        rendered = f"{number:.6e}"
    return f"{rendered} {latex_escape(unit)}".strip()


def rows(items: list[tuple[str, str, str]]) -> str:
    return "\n".join(
        rf"{latex_escape(left)} & {latex_escape(mid)} & {right} \\"
        for left, mid, right in items
    )


def comparison_rows() -> str:
    items = [
        (
            "BeamBeam field dumps",
            "sandbox/BeamBeam-static-1V.in",
            "ASCII rho/phi/E dumps in data/sandbox are compared against the analytic symmetric Gaussian manufactured solution.",
        ),
        (
            "Gamma-gamma e-/e+ tracking",
            "sandbox/track-e-p/gamma_gamma_pairs-2.in",
            "FROMFILE pair kinematics and final witness-container .stat summaries are checked for the primary beam, electrons, and positrons.",
        ),
        (
            "OPALX-IMPACT drift",
            "sandbox/OPALX-IMPACT/drift-1.in",
            "A simple proton drift with space charge is compared with digitized IMPACT reference data in extracted_graph_data.csv.",
        ),
    ]
    return "\n".join(
        rf"\textbf{{{latex_escape(a)}}} & {{\footnotesize\path{{{b}}}}} & {latex_escape(c)} \\"
        for a, b, c in items
    )


def selected_tables(metrics: dict[str, object]) -> str:
    field_items = [
        ("field snapshot", "diagnostic kind", fmt(metrics, "fields.snapshot_kind")),
        ("field centers", "number of source centers", fmt(metrics, "fields.n_centers")),
        ("charge integral", "deposited charge", fmt(metrics, "fields.charge_c", "C")),
        ("rho relative L2", "density error", fmt(metrics, "fields.rho.rel_l2")),
        ("phi relative L2", "potential error", fmt(metrics, "fields.phi.rel_l2")),
        ("Ex relative L2", "electric field error", fmt(metrics, "fields.Ex.rel_l2")),
        ("Ey relative L2", "electric field error", fmt(metrics, "fields.Ey.rel_l2")),
        ("Ez relative L2", "electric field error", fmt(metrics, "fields.Ez.rel_l2")),
    ]
    impact_items = [
        ("OPALX rows", "stat samples", fmt(metrics, "impact.opalx_rows")),
        ("IMPACT rows", "digitized samples", fmt(metrics, "impact.digitized_rows")),
        ("final x", "OPALX final centroid", fmt(metrics, "impact.final.OPALX_x_m", "m")),
        ("final epsx", "OPALX final emittance", fmt(metrics, "impact.final.OPALX_epsx_mm_mr", "mm mrad")),
        ("RMS x difference", "OPALX - IMPACT", fmt(metrics, "impact.rms_x_diff_m", "m")),
        ("max x difference", "OPALX - IMPACT", fmt(metrics, "impact.max_abs_x_diff_m", "m")),
        ("RMS epsx difference", "OPALX - IMPACT", fmt(metrics, "impact.rms_epsx_diff_mm_mr", "mm mrad")),
        ("max epsx difference", "OPALX - IMPACT", fmt(metrics, "impact.max_abs_epsx_diff_mm_mr", "mm mrad")),
    ]
    pair_items = [
        ("electron count", "FROMFILE particles", fmt(metrics, "pairs.fromfile.electron.N")),
        ("positron count", "FROMFILE particles", fmt(metrics, "pairs.fromfile.positron.N")),
        ("electron mean K", "input pair kinematics", fmt(metrics, "pairs.fromfile.electron.<K> keV", "keV")),
        ("positron mean K", "input pair kinematics", fmt(metrics, "pairs.fromfile.positron.<K> keV", "keV")),
        ("median opening angle", "pair summary", fmt(metrics, "pairs.pair_summary.median opening angle deg", "deg")),
        ("opening > 150 deg", "pair summary fraction", fmt(metrics, "pairs.pair_summary.fraction opening > 150 deg")),
        ("electron final energy", "gamma_gamma_pairs-2 c1", fmt(metrics, "pairs.track_stat.gamma_gamma_pairs-2.c1.final_energy_MeV", "MeV")),
        ("positron final energy", "gamma_gamma_pairs-2 c2", fmt(metrics, "pairs.track_stat.gamma_gamma_pairs-2.c2.final_energy_MeV", "MeV")),
        ("electron final rms x", "gamma_gamma_pairs-2 c1", fmt(metrics, "pairs.track_stat.gamma_gamma_pairs-2.c1.final_rms_x_m", "m")),
        ("positron final rms x", "gamma_gamma_pairs-2 c2", fmt(metrics, "pairs.track_stat.gamma_gamma_pairs-2.c2.final_rms_x_m", "m")),
    ]
    template = r"""
\subsection*{BeamBeam Field Manufactured Solution}
\begin{tabularx}{\textwidth}{@{}p{0.25\textwidth}p{0.35\textwidth}X@{}}
\toprule
Metric & Meaning & Current value \\
\midrule
%s
\bottomrule
\end{tabularx}

\subsection*{OPALX-IMPACT Drift}
\begin{tabularx}{\textwidth}{@{}p{0.25\textwidth}p{0.35\textwidth}X@{}}
\toprule
Metric & Meaning & Current value \\
\midrule
%s
\bottomrule
\end{tabularx}

\begin{figure}[h]
\centering
\includegraphics[width=0.82\textwidth]{opalx_impact_drift_comparison.png}
\caption{Digitized IMPACT reference data overlaid with the current OPALX drift result from \texttt{sandbox/OPALX-IMPACT/drift-1.stat}.}
\end{figure}

\subsection*{Gamma-Gamma Pair Tracking}
\begin{tabularx}{\textwidth}{@{}p{0.25\textwidth}p{0.35\textwidth}X@{}}
\toprule
Metric & Meaning & Current value \\
\midrule
%s
\bottomrule
\end{tabularx}
"""
    return template % (rows(field_items), rows(impact_items), rows(pair_items))


def collision_timeline_figure() -> str:
    return r"""
\begin{figure}[h]
\centering
\begin{tikzpicture}[
    x=1cm,
    y=-0.065cm,
    >=Latex,
    font=\footnotesize,
    bunch/.style={draw=black, line width=0.45pt, minimum width=0.58cm, minimum height=0.22cm},
    event/.style={anchor=west, align=left, text width=5.6cm},
    window/.style={draw=black, line width=0.45pt, minimum width=2.5cm, minimum height=0.62cm},
]
\definecolor{sourcegreen}{HTML}{1B9E77}
\definecolor{pairmagenta}{HTML}{D33682}
\definecolor{paircyan}{HTML}{2AA6B8}

\draw[->, line width=0.6pt] (0,15) -- (0,212);
\node[anchor=south] at (0,12) {time [ps]};
\foreach \t/\y in {22/22,45/45,100/100,111/136,116/169,121/202} {
    \draw[line width=0.45pt] (-0.12,\y) -- (0.12,\y);
    \node[left] at (-0.18,\y) {\t};
}

\node[event] at (8.2,22) {approach: source bunches still separated};
\node[event] at (8.2,45) {\texttt{BB-active=TRUE}};
\node[event] at (8.2,100) {witness pair born at IP: \texttt{c1,c2:n=1297}};
\node[event] at (8.2,136) {\texttt{source\_bunches\_overlap=TRUE}};
\node[event] at (8.2,169) {source bunches passed IP};
\node[event] at (8.2,202) {\texttt{source\_bunches\_overlap=FALSE}\\
\texttt{retired\_bunches=1, source\_active=FALSE}};

\foreach \y/\xb/\xg in {22/3.55/6.45,45/4.35/5.65,100/4.78/5.22,136/5.00/5.00,169/5.35/4.65} {
    \node[window] at (5,\y) {};
    \draw[line width=0.45pt] (5,\y-6.0) -- (5,\y+6.0);
    \node[bunch, fill=black] at (\xb,\y) {};
    \node[bunch, draw=sourcegreen, fill=sourcegreen!18, line width=0.85pt] at (\xg,\y) {};
}
\node[window] at (5,202) {};
\draw[line width=0.45pt] (5,196.0) -- (5,208.0);
\node[bunch, draw=gray, dashed] (retiredPrimaryA) at (4.65,202) {};
\node[bunch, draw=gray, dashed] (retiredPrimaryB) at (5.35,202) {};
\foreach \nodeName in {retiredPrimaryA,retiredPrimaryB} {
    \draw[gray, line width=0.55pt]
        ([xshift=-5pt,yshift=-3pt]\nodeName.center) --
        ([xshift=5pt,yshift=3pt]\nodeName.center);
    \draw[gray, line width=0.55pt]
        ([xshift=-5pt,yshift=3pt]\nodeName.center) --
        ([xshift=5pt,yshift=-3pt]\nodeName.center);
}

\foreach \t/\h in {100/8,136/10,169/11,202/12} {
    \path[draw=black, fill=paircyan!55, line width=0.45pt]
        (5.00,\t) -- ([xshift=-4pt,yshift=-\h pt]5.00,\t)
        -- ([xshift=4pt,yshift=-\h pt]5.00,\t) -- cycle;
    \path[draw=black, fill=pairmagenta!65, line width=0.45pt]
        (5.00,\t) -- ([xshift=-4pt,yshift=\h pt]5.00,\t)
        -- ([xshift=4pt,yshift=\h pt]5.00,\t) -- cycle;
}

\node[anchor=south] at (5.0,11) {BeamBeam window};
\node[anchor=north] at (5.24,16) {IP};
\draw[->, line width=0.45pt] (3.45,16) -- (4.05,16);
\draw[->, line width=0.45pt] (6.55,16) -- (5.95,16);
\node[anchor=west, text width=8.2cm, font=\scriptsize] at (0.7,217) {Draft event timeline, not to scale. The rows are meant to track grepable \texttt{BB-DIAG} transitions from the OPALX stdout; the final row shows the configured \texttt{RETIRE\_TIME} deleting the primary source while the witness pair remains active.};
\end{tikzpicture}
\caption{Draft timeline for the copied BeamBeam source bunches and the injected low-energy witness pair.}
\end{figure}
"""


def large_cylinder_timeline_figure() -> str:
    return r"""
\begin{figure}[h]
\centering
\begin{tikzpicture}[x=1cm, y=1cm, font=\footnotesize]
\definecolor{sourcegreen}{HTML}{1B9E77}
\definecolor{pairmagenta}{HTML}{D33682}
\definecolor{paircyan}{HTML}{2AA6B8}
\draw[line width=0.7pt] (0,0.55) -- (8.6,0.55);
\draw[line width=0.7pt] (0,-0.55) -- (8.6,-0.55);
\draw[line width=0.7pt] (0,0) ellipse (0.32 and 0.55);
\draw[line width=0.7pt] (8.6,0) ellipse (0.32 and 0.55);
\draw[line width=0.5pt] (4.3,-0.78) -- (4.3,0.78);
\node[anchor=north] at (4.3,-0.82) {IP, $s=0.17$ m};
\node[anchor=south] at (4.3,0.68) {BeamBeam cylinder: length 32 cm, radius 15 cm};
\node[draw=black, fill=paircyan!55, minimum width=0.20cm, minimum height=0.40cm] at (4.18,0.0) {};
\node[draw=black, fill=pairmagenta!65, minimum width=0.20cm, minimum height=0.40cm] at (4.42,0.0) {};
\end{tikzpicture}

\vspace{0.6em}
\begin{tabularx}{0.86\textwidth}{@{}p{0.38\textwidth}p{0.18\textwidth}X@{}}
\toprule
Event & Time & Note \\
\midrule
BeamBeam window active & observed & \texttt{BB-active=TRUE} is emitted before witness injection. \\
Witness pair injected & 550.3 ps & Preserves the 16.747 ps pre-IP offset from the proof-of-principle setup. \\
Primary reaches IP & 567.1 ps & $s_\mathrm{IP}=0.17$ m for the 32 cm BeamBeam window. \\
Source retired & 571.3 ps & Preserves the 4.253 ps post-IP retirement offset. \\
Tracking ends & 1600 ps & Long enough for the first 15 cm radial crossings. \\
\bottomrule
\end{tabularx}
\caption{Timing and geometry summary for \texttt{gamma\_gamma\_pairs-large-cylinder.in}.}
\end{figure}
"""


def baseline_provenance_block(metadata: dict[str, object]) -> str:
    git = metadata.get("git", {}) if isinstance(metadata, dict) else {}
    if not isinstance(git, dict) or not git:
        return "No baseline git metadata is available. Regenerate the accepted baseline with \\texttt{--update-baseline}."

    status = str(git.get("status_short", "")).strip()
    if not status:
        status_text = "clean"
    else:
        status_text = r"\\".join(latex_escape(line) for line in status.splitlines())
    commit = str(git.get("commit", "unknown"))
    rows_text = rows(
        [
            ("branch", "baseline source", latex_escape(git.get("branch", "unknown"))),
            ("commit", "baseline source", rf"\texttt{{{latex_escape(commit)}}}"),
            ("dirty", "baseline source", latex_escape(str(git.get("dirty", "unknown")))),
            ("status", "git status --short when accepted", rf"\begin{{minipage}}[t]{{0.38\textwidth}}\raggedright\ttfamily\footnotesize {status_text}\end{{minipage}}"),
        ]
    )
    return rf"""\begin{{tabularx}}{{\textwidth}}{{@{{}}p{{0.22\textwidth}}p{{0.28\textwidth}}X@{{}}}}
\toprule
Item & Meaning & Value \\
\midrule
{rows_text}
\bottomrule
\end{{tabularx}}"""


def document(
        metrics: dict[str, object],
        csv_path: Path,
        baseline_path: Path,
        metadata: dict[str, object],
        loss_plot_path: Path,
        large_cylinder_plot_path: Path,
        large_cylinder_charge_compare_plot_path: Path,
) -> str:
    generated = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    return rf"""\documentclass[10pt]{{article}}
\usepackage[margin=0.75in]{{geometry}}
\usepackage{{array}}
\usepackage{{booktabs}}
\usepackage{{graphicx}}
\usepackage[colorlinks=true,linkcolor=black,urlcolor=blue]{{hyperref}}
\usepackage{{longtable}}
\usepackage{{tabularx}}
\usepackage{{tikz}}
\usepackage[table]{{xcolor}}
\usetikzlibrary{{arrows.meta}}
\definecolor{{opalxblue}}{{HTML}}{{214B78}}
\definecolor{{opalxlight}}{{HTML}}{{EEF4FA}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.45em}}
\renewcommand{{\arraystretch}}{{1.18}}

\title{{\textbf{{OPALX BeamBeam Sandbox Regression Overview}}}}
\author{{Accepted baseline: \texttt{{{latex_escape(str(baseline_path.relative_to(ROOT)))}}}\\Latest metrics: \texttt{{{latex_escape(str(csv_path.relative_to(ROOT)))}}}}}
\date{{{latex_escape(generated)}}}

\begin{{document}}
\maketitle

\begin{{center}}
\colorbox{{opalxlight}}{{\parbox{{0.94\textwidth}}{{This overview summarizes the sandbox regression checks.  The accepted baseline values are stored in \texttt{{{latex_escape(str(baseline_path.relative_to(ROOT)))}}}; this is the comparison reference.  The latest run is stored separately in \texttt{{{latex_escape(str(csv_path.relative_to(ROOT)))}}}; it is regenerated by the harness and is not the source of truth.  The accepted baseline records its repository state in \texttt{{\_metadata.git}}, shown below.}}}}
\end{{center}}

\section*{{What Is Compared}}
\rowcolors{{2}}{{opalxlight}}{{white}}
\begin{{tabularx}}{{\textwidth}}{{@{{}}p{{0.23\textwidth}}p{{0.28\textwidth}}X@{{}}}}
\toprule
Check & Input / source & Comparison target \\
\midrule
{comparison_rows()}
\bottomrule
\end{{tabularx}}
\rowcolors{{2}}{{white}}{{white}}

\section*{{Accepted Baseline Provenance}}
The baseline values were accepted from the git state below.  If \texttt{{dirty=True}}, the commit alone is not sufficient to reproduce the baseline; the listed \texttt{{git status --short}} entries are part of the provenance.

{baseline_provenance_block(metadata)}

\section*{{Physics Context From The Note}}
The existing note \texttt{{sandbox/note/boosted\_gaussian\_witness.tex}} defines the manufactured Gaussian witness model.  The overview reuses its generated figures as visual context; the note files themselves are not part of the regression baseline.

\begin{{figure}}[h]
\centering
\includegraphics[width=0.48\textwidth]{{../note/figs/boosted_gaussian_witness_physical_field_t0.png}}
\hfill
\includegraphics[width=0.48\textwidth]{{../note/figs/boosted_gaussian_witness_physical_paths.png}}
\caption{{Field and trajectory context generated by the BeamBeam witness note.}}
\end{{figure}}

\section*{{Current Regression Metrics}}
{selected_tables(metrics)}

\section*{{BeamBeam Diagnostic Timeline}}
{collision_timeline_figure()}

\section*{{Cylinder Loss Timeline}}
\begin{{figure}}[h]
\centering
\includegraphics[width=0.88\textwidth]{{{loss_plot_path.name}}}
\caption{{Cylinder-edge crossing histogram for the low-energy witness pair after source retirement.  Bars show the first recorded H5 particle position with radius at least 1 cm for containers \texttt{{c1}} and \texttt{{c2}}; the lower panel sketches the cylindrical BeamBeam aperture.  In this run OPALX records the particle trajectories but does not delete the witness particles at the BeamBeam cylinder boundary.}}
\end{{figure}}

\clearpage
\section*{{Large Cylinder Setup}}
{large_cylinder_timeline_figure()}

\begin{{figure}}[h]
\centering
\includegraphics[width=0.88\textwidth]{{{large_cylinder_plot_path.name}}}
\caption{{Cylinder-edge crossing histogram for the large setup with 15 cm radius and 32 cm BeamBeam length.  Bars show first recorded H5 particle positions with radius at least 15 cm for containers \texttt{{c1}} and \texttt{{c2}}.  The lower panel marks the nominal 32 cm BeamBeam cylinder; crossings outside that longitudinal span are shown rather than clipped.}}
\end{{figure}}

\begin{{figure}}[h]
\centering
\includegraphics[width=0.88\textwidth]{{{large_cylinder_charge_compare_plot_path.name}}}
\caption{{Comparison of the nominal large-cylinder crossing histograms with a rerun where the primary bunch charge is scaled by $10^{{-5}}$.  Solid colored curves show the nominal primary charge, dashed black curves show the reduced-charge run, and grey bars show the reduced-minus-nominal bin difference.  Both runs give 701 electron and 692 positron first crossings, with zero bin-count difference in this histogram; matching particle IDs shift by at most about 141--142 $\mu$m in first-crossing $z$.}}
\end{{figure}}

\section*{{How To Reproduce}}
\begin{{verbatim}}
source .venv-h6/bin/activate
python sandbox/python/run_sandbox_regressions.py \
  --run-opalx \
  --opalx-exe build_openmp/src/opalx
cd sandbox/track-e-p
../../build_openmp/src/opalx gamma_gamma_pairs-large-cylinder.in
../../build_openmp/src/opalx gamma_gamma_pairs-large-cylinder-q1em5.in
cd ../..
python sandbox/python/make_sandbox_regression_overview.py
\end{{verbatim}}

\end{{document}}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--impact-plot", type=Path, default=DEFAULT_IMPACT_PLOT)
    parser.add_argument("--loss-plot", type=Path, default=DEFAULT_LOSS_PLOT)
    parser.add_argument("--large-cylinder-plot", type=Path, default=DEFAULT_LARGE_CYLINDER_PLOT)
    parser.add_argument("--large-cylinder-charge-compare-plot", type=Path, default=DEFAULT_LARGE_CYLINDER_CHARGE_COMPARE_PLOT)
    parser.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--no-pdf", action="store_true", help="Only write the LaTeX source.")
    parser.add_argument("--skip-impact-plot", action="store_true", help="Reuse the existing OPALX-IMPACT plot.")
    parser.add_argument("--skip-loss-plot", action="store_true", help="Reuse the existing cylinder loss plot.")
    parser.add_argument("--skip-large-cylinder-plot", action="store_true", help="Reuse the existing large-cylinder crossing plot.")
    parser.add_argument("--skip-large-cylinder-charge-compare-plot", action="store_true", help="Reuse the existing large-cylinder charge-comparison plot.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = metric_map(args.csv)
    metadata = baseline_metadata(args.baseline)
    if not args.skip_impact_plot:
        make_impact_plot(args.impact_plot)
    if not args.skip_loss_plot:
        make_cylinder_loss_plot(args.loss_plot)
    if not args.skip_large_cylinder_plot:
        make_large_cylinder_crossing_plot(args.large_cylinder_plot)
    if not args.skip_large_cylinder_charge_compare_plot:
        make_large_cylinder_charge_compare_plot(args.large_cylinder_charge_compare_plot)
    args.tex.write_text(
            document(
                    metrics, args.csv, args.baseline, metadata, args.loss_plot,
                    args.large_cylinder_plot, args.large_cylinder_charge_compare_plot),
            encoding="utf-8")
    print(f"wrote {args.tex}")

    if not args.no_pdf:
        subprocess.run(
            ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", args.tex.name],
            cwd=args.tex.parent,
            check=True,
        )
        print(f"wrote {args.tex.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
