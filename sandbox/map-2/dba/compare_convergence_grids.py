#!/usr/bin/env python3
"""Pair preserved before/after DBA studies without changing either study's data.

Runs no simulation. Reports all matched cases, full-map-selected optima, and
publication-quality comparison figures. A gain greater than one means less error.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from convergence_dt import np, pd, plt
from convergence_grid_plots import save
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import SymmetricalLogLocator

KEYS = ["dt_s", "epsilon", "richardson_levels", "ranks"]
METRICS = ["scaled_max_error", "abs_R16", "abs_R26"]
LABELS = [r"Full-map max-entry error", r"$|R_{16}|$ [m]", r"$|R_{26}|$"]
STYLES = [("before", "Before", "#757575", "o", "--"),
          ("after", "After bookkeeping", "#0072B2", "s", "-")]


def load_paired(before: Path, after: Path) -> pd.DataFrame:
    manifests = [json.loads((root / "provenance.json").read_text()) for root in (before, after)]
    for field in ("input_template_sha256", "distribution_sha256"):
        if manifests[0]["identity"][field] != manifests[1]["identity"][field]:
            raise ValueError(f"incompatible physical inputs: {field}")
    if not np.array_equal(np.loadtxt(before / "analytic-map.txt"), np.loadtxt(after / "analytic-map.txt")):
        raise ValueError("analytic references differ")
    frames = [pd.read_csv(root / "results.csv", float_precision="round_trip") for root in (before, after)]
    primary = [frame[frame.ranks == 1].copy() for frame in frames]
    if set(map(tuple, primary[0][KEYS].to_numpy())) != set(map(tuple, primary[1][KEYS].to_numpy())):
        raise ValueError("primary grids differ")
    paired = primary[0].merge(primary[1], on=KEYS, how="inner", suffixes=("_before", "_after"),
                             validate="one_to_one").copy()
    if not (paired.status_before == paired.status_after).all():
        raise ValueError("acceptance/rejection behavior changed")
    for _, row in paired.iterrows():
        old = json.loads((before / row.path_before / "result.json").read_text())
        new = json.loads((after / row.path_after / "result.json").read_text())
        if old["artifacts_sha256"]["map-2-dba.in"] != new["artifacts_sha256"]["map-2-dba.in"]:
            raise ValueError(f"configured inputs differ: {row.path_before}")
    for metric in METRICS:
        paired[f"{metric}_gain"] = paired[f"{metric}_before"] / paired[f"{metric}_after"]
    return paired


def plot_pair(axis, data, x, metric, signed=False):
    data = data.sort_values(x)
    for suffix, label, color, marker, style in STYLES:
        axis.plot(data[x], data[f"{metric}_{suffix}"], label=label, color=color,
                  marker=marker, linestyle=style, ms=4.0, lw=1.4)
    axis.set_xscale("log")
    axis.invert_xaxis()
    if signed:
        axis.set_yscale("symlog", linthresh=1e-9)
        locator = SymmetricalLogLocator(linthresh=1e-9, base=10)
        locator.set_params(numticks=7)
        axis.yaxis.set_major_locator(locator)
        axis.axhline(0, color="0.35", lw=.7, ls=":")
    else:
        axis.set_yscale("log")
    axis.set_xlabel(r"Time step $\Delta t$ [s]" if x == "dt_s" else r"Starting $\epsilon$")
    axis.grid(which="major", alpha=.25, lw=.6)


def finish(figure, output, name, title, footer):
    figure.suptitle(title, fontsize=14)
    handles = [Line2D([], [], color=color, marker=marker, ls=style, label=label)
               for _, label, color, marker, style in STYLES]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .953),
                  ncol=2, frameon=False)
    figure.get_layout_engine().set(rect=(0, .075, 1, .81))
    figure.text(.5, .022, footer, ha="center", fontsize=9)
    save(figure, output, name)


def make_comparison_plots(valid, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 11,
                         "xtick.labelsize": 9, "ytick.labelsize": 9, "svg.fonttype": "none",
                         "axes.spines.top": False, "axes.spines.right": False})
    l1 = valid[valid.richardson_levels == 1]
    cuts = [(l1[l1.epsilon == .003], "dt_s", r"L=1, $\epsilon=3\times10^{-3}$"),
            (l1[l1.dt_s == 1e-13], "epsilon", r"L=1, $\Delta t=10^{-13}$ s")]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.5), layout="constrained")
    for row, (data, x, title) in enumerate(cuts):
        for column, (metric, label) in enumerate(zip(METRICS, LABELS)):
            plot_pair(axes[row, column], data, x, metric)
            axes[row, column].set_ylabel(label)
            axes[row, column].set_title(title, fontsize=11)
    finish(fig, output, "before-after", "DBA: same settings, before and after bookkeeping",
           r"Boris unchanged; $\epsilon$ in m for $(x,y,\zeta)$, dimensionless otherwise. Analytic $R_{16}=R_{26}=0$.")

    fig, axes = plt.subplots(5, 3, figsize=(11.5, 14.0), layout="constrained")
    for level in range(5):
        data = valid[(valid.richardson_levels == level) & (valid.dt_s == 1e-13)]
        for column, (metric, label) in enumerate(zip(METRICS, LABELS)):
            plot_pair(axes[level, column], data, "epsilon", metric)
            axes[level, column].set_ylabel(label)
            axes[level, column].set_title(f"Richardson L={level}", fontsize=11)
    finish(fig, output, "all-levels-vs-epsilon", r"DBA: every refinement level at $\Delta t=10^{-13}$ s",
           "All paired amplitudes shown; the integrator and finite-difference formulas are unchanged.")

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.5), layout="constrained")
    for row, (data, x, title) in enumerate(cuts):
        for column, metric in enumerate(["R16", "R26"]):
            plot_pair(axes[row, column], data, x, metric, signed=True)
            axes[row, column].set_ylabel(r"$R_{16}$ [m]" if column == 0 else r"$R_{26}$")
            axes[row, column].set_title(title, fontsize=11)
    finish(fig, output, "signed-before-after", "DBA: signed dispersion entries at matched settings",
           r"Symmetric-log vertical axes, linear within $\pm10^{-9}$; analytic values and dotted line = zero.")

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.3), layout="constrained")
    dts = sorted(valid.dt_s.unique(), reverse=True)
    epsilons = sorted(valid.epsilon.unique(), reverse=True)
    extent = max(1.0, np.ceil(np.abs(np.log10(valid.scaled_max_error_gain)).max()))
    for level, axis in zip(range(5), axes.flat):
        data = valid[valid.richardson_levels == level]
        grid = data.pivot(index="dt_s", columns="epsilon", values="scaled_max_error_gain").reindex(
            index=dts, columns=epsilons)
        image = axis.imshow(np.log10(grid.to_numpy()), aspect="auto", cmap="RdBu",
                            norm=TwoSlopeNorm(vmin=-extent, vcenter=0, vmax=extent))
        axis.set_xticks(range(len(epsilons)), [f"{value:.0e}" for value in epsilons], rotation=60)
        axis.set_yticks(range(len(dts)), [f"{value:.0e}" for value in dts])
        axis.tick_params(labelsize=8)
        axis.set_xlabel(r"Starting $\epsilon$")
        axis.set_ylabel(r"$\Delta t$ [s]")
        axis.set_title(f"Richardson L={level}")
    axes[1, 2].set_axis_off()
    fig.colorbar(image, ax=list(axes.flat), shrink=.8, pad=.025,
                 label=r"$\log_{10}(E_{\rm before}/E_{\rm after})$: positive = improvement")
    fig.suptitle("DBA: full-map accuracy gain across all 330 valid paired cases", fontsize=14)
    save(fig, output, "paired-full-map-gains")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in (args.before.resolve(), args.after.resolve()):
        parser.error("use a separate comparison directory")
    paired = load_paired(args.before, args.after)
    args.output.mkdir(parents=True, exist_ok=True)
    paired.to_csv(args.output / "paired-results.csv", index=False, float_format="%.12e")
    valid = paired[paired.status_before == "OK"]
    best = []
    for level, group in valid.groupby("richardson_levels"):
        item = {"richardson_levels": level}
        for suffix in ["before", "after"]:
            row = group.loc[group[f"scaled_max_error_{suffix}"].idxmin()]
            item.update({f"{key}_{suffix}": row[key] for key in ["dt_s", "epsilon"]})
            item.update({f"{key}_{suffix}": row[f"{key}_{suffix}"] for key in
                         ["scaled_max_error", "R16", "R26", "relative_frobenius_error", "worst_entry"]})
        item["best_error_gain"] = item["scaled_max_error_before"] / item["scaled_max_error_after"]
        best.append(item)
    best = pd.DataFrame(best)
    best.to_csv(args.output / "best-before-after.csv", index=False, float_format="%.12e")
    valid.groupby(["richardson_levels", "dt_s"])[[f"{m}_gain" for m in METRICS]].median().to_csv(
        args.output / "median-gains.csv", float_format="%.12e")
    summary = {"primary_cases": len(paired), "valid_pairs": len(valid),
               "improved_pairs": int((valid.scaled_max_error_gain > 1).sum()),
               "worsened_pairs": int((valid.scaled_max_error_gain < 1).sum()),
               "median_full_map_gain": float(valid.scaled_max_error_gain.median()),
               "gain_definition": "before error / after error; >1 improves; not a significance test",
               "best_before": float(valid.scaled_max_error_before.min()),
               "best_after": float(valid.scaled_max_error_after.min()),
               "input_checks": "same template, distribution, analytic map, all per-case input hashes",
               "sources": {str(root.resolve()): hashlib.sha256((root / "results.csv").read_bytes()).hexdigest()
                           for root in (args.before, args.after)},
               "comparison_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_comparison_plots(valid, args.output)
    print(best.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
