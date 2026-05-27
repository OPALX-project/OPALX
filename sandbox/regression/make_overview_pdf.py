#!/usr/bin/env python3
"""Build a compact PDF overview of the sandbox regression comparisons."""

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
DEFAULT_CSV = HERE / "current_metrics.csv"
DEFAULT_TEX = HERE / "sandbox_regression_overview.tex"
DEFAULT_BASELINE = HERE / "sandbox_regression_baseline.json"
DEFAULT_IMPACT_PLOT = HERE / "opalx_impact_drift_comparison.png"


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
\usepackage[table]{{xcolor}}
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

\section*{{How To Reproduce}}
\begin{{verbatim}}
source .venv-h6/bin/activate
python sandbox/regression/run_sandbox_regressions.py \
  --run-opalx \
  --opalx-exe build_openmp/src/opalx
python sandbox/regression/make_overview_pdf.py
\end{{verbatim}}

\end{{document}}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--impact-plot", type=Path, default=DEFAULT_IMPACT_PLOT)
    parser.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--no-pdf", action="store_true", help="Only write the LaTeX source.")
    parser.add_argument("--skip-impact-plot", action="store_true", help="Reuse the existing OPALX-IMPACT plot.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = metric_map(args.csv)
    metadata = baseline_metadata(args.baseline)
    if not args.skip_impact_plot:
        make_impact_plot(args.impact_plot)
    args.tex.write_text(document(metrics, args.csv, args.baseline, metadata), encoding="utf-8")
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
