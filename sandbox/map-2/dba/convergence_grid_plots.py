#!/usr/bin/env python3
"""Static publication-quality figures for the recorded DBA convergence grid."""

from __future__ import annotations

from pathlib import Path

from convergence_dt import np, pd, plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import SymmetricalLogLocator

COLORS = ["#333333", "#0072B2", "#D55E00", "#009E73", "#CC79A7"]
MARKERS = ["o", "s", "^", "D", "v"]
METRICS = ["scaled_max_error", "abs_R16", "abs_R26"]
LABELS = [r"Full-map error $E_\infty$", r"$|R_{16}|$ [m]", r"$|R_{26}|$"]
LEVEL_LABELS = ["Centered (L=0)", "Richardson L=1", "Richardson L=2",
                "Richardson L=3", "Richardson L=4"]


def scientific(value: float) -> str:
    exponent = int(np.floor(np.log10(abs(value))))
    coefficient = value / 10**exponent
    return (rf"10^{{{exponent}}}" if np.isclose(coefficient, 1)
            else rf"{coefficient:g}\times10^{{{exponent}}}")


def save(figure, output: Path, name: str) -> None:
    figure.savefig(output / f"{name}.png", dpi=220, facecolor="white")
    figure.savefig(output / f"{name}.svg", facecolor="white")
    plt.close(figure)


def legend(figure) -> None:
    handles = [Line2D([], [], color=color, marker=marker, lw=1.5, ms=5, label=label)
               for color, marker, label in zip(COLORS, MARKERS, LEVEL_LABELS)]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.943),
                  ncol=5, frameon=False, fontsize=9)


def lines(axis, data: pd.DataFrame, x: str, metric: str, signed: bool = False) -> None:
    for level, group in data.groupby("richardson_levels"):
        level = int(level)
        group = group.sort_values(x)
        y = group[metric] if signed else group[metric].where(group[metric] > 0)
        axis.plot(group[x], y, color=COLORS[level], marker=MARKERS[level],
                  lw=1.4, ms=4, label=LEVEL_LABELS[level])
    axis.set_xscale("log")
    if signed:
        axis.set_yscale("symlog", linthresh=1e-9, linscale=1.0)
        locator = SymmetricalLogLocator(linthresh=1e-9, base=10)
        locator.set_params(numticks=9)
        axis.yaxis.set_major_locator(locator)
    else:
        axis.set_yscale("log")
    axis.invert_xaxis()
    axis.grid(True, which="major", lw=0.6, alpha=0.24)
    axis.grid(True, which="minor", lw=0.4, alpha=0.08)
    if signed:
        axis.axhline(0, color="0.35", ls=":", lw=0.8)


def cuts(primary: pd.DataFrame, output: Path, x: str, fixed: str,
         values: list[float], name: str) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(11.4, 7.3), layout="constrained")
    figure.set_constrained_layout_pads(h_pad=0.07, w_pad=0.07, hspace=0.08, wspace=0.03)
    for row, value in enumerate(values):
        data = primary[primary[fixed] == value]
        for column, (metric, ylabel) in enumerate(zip(METRICS, LABELS)):
            axis = axes[row, column]
            lines(axis, data, x, metric)
            axis.set_ylabel(ylabel)
            axis.set_xlabel(r"Time step $\Delta t$ [s]" if x == "dt_s" else r"Starting amplitude $\epsilon$")
            title = (rf"$\epsilon={scientific(value)}$" if fixed == "epsilon"
                     else rf"$\Delta t={scientific(value)}$ s")
            axis.set_title(title, fontsize=11)
    if x == "dt_s":
        title = "DBA: time-step convergence at two fixed starting amplitudes"
        footer = r"Boris/LF2 with boundary splitting; $\Delta t=10^{-10}$ s rejected by the existing QACH guard."
    else:
        title = "DBA: perturbation convergence at two fixed time steps"
        footer = r"$\epsilon$ is in m for $(x,y,\zeta)$ and dimensionless for $(x',y',\delta)$; analytic $R_{16}=R_{26}=0$."
    figure.suptitle(title, fontsize=14)
    # Reserve a clean legend strip separate from the constrained-layout axes.
    figure.get_layout_engine().set(rect=(0, 0.065, 1, 0.835))
    legend(figure)
    figure.text(0.5, 0.02, footer, ha="center", fontsize=9)
    save(figure, output, name)


def heatmaps(frame: pd.DataFrame, output: Path) -> None:
    primary = frame[frame.ranks == 1]
    dts = sorted(primary.dt_s.unique(), reverse=True)
    epsilons = sorted(primary.epsilon.unique(), reverse=True)
    levels = sorted(primary.richardson_levels.unique())
    valid = primary[primary.status == "OK"]
    values = valid[METRICS].to_numpy().flatten()
    values = values[np.isfinite(values) & (values > 0)]
    # Same numeric color scale for every panel; the column titles state units.
    norm = LogNorm(10 ** np.floor(np.log10(values.min())), 10 ** np.ceil(np.log10(values.max())))
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("#d9d9d9")
    figure, axes = plt.subplots(len(levels), 3, figsize=(11.8, 14.2), layout="constrained")
    for row, level in enumerate(levels):
        subset = primary[primary.richardson_levels == level]
        for column, metric in enumerate(METRICS):
            grid = subset.pivot(index="dt_s", columns="epsilon", values=metric).reindex(
                index=dts, columns=epsilons)
            image = axes[row, column].imshow(np.ma.masked_invalid(grid.to_numpy()),
                                            norm=norm, cmap=cmap, aspect="auto")
            axes[row, column].set_yticks(range(len(dts)), [f"{dt:.0e}" for dt in dts])
            axes[row, column].set_xticks(range(len(epsilons)),
                                        [f"{eps:.0e}" for eps in epsilons], rotation=45)
            axes[row, column].tick_params(labelsize=8)
            axes[row, column].set_xlabel(r"Starting $\epsilon$", fontsize=9)
            axes[row, column].set_ylabel(r"$\Delta t$ [s]", fontsize=9)
            axes[row, column].set_title(f"L={level}: " + LABELS[column], fontsize=11)
            if metric in subset and subset[metric].notna().any():
                best = subset.loc[subset[metric].idxmin()]
                axes[row, column].plot(epsilons.index(best.epsilon), dts.index(best.dt_s),
                                       "*", ms=12, color="white", markeredgecolor="black", mew=0.7)
    figure.colorbar(image, ax=axes, shrink=0.65, pad=0.02, label="Error / residual magnitude")
    figure.suptitle("DBA: complete DT × amplitude × Richardson grid\n"
                    "Gray = rejected DT; star = smallest residual in that panel", fontsize=14)
    save(figure, output, "complete-grid")


def best_vs_dt(primary: pd.DataFrame, output: Path) -> None:
    best = primary.loc[primary.groupby(["richardson_levels", "dt_s"]).scaled_max_error.idxmin()]
    figure, axes = plt.subplots(1, 3, figsize=(11.4, 4.5), layout="constrained")
    for axis, metric, ylabel in zip(axes, METRICS, LABELS):
        lines(axis, best, "dt_s", metric)
        axis.set_xlabel(r"Time step $\Delta t$ [s]")
        axis.set_ylabel(ylabel)
    figure.suptitle("DBA: best tested starting amplitude at each time step", fontsize=14)
    figure.get_layout_engine().set(rect=(0, 0.1, 1, 0.76))
    legend(figure)
    figure.text(0.5, 0.025, r"Each point selects $\epsilon$ by the full-map error; "
                r"$R_{16}$ and $R_{26}$ use that same map, not separate optima.", ha="center", fontsize=9)
    save(figure, output, "best-vs-dt")


def signed_entries(primary: pd.DataFrame, output: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(9.8, 7.6), layout="constrained")
    slices = [(primary[primary.epsilon == 1e-3], "dt_s", r"$\epsilon=10^{-3}$"),
              (primary[primary.dt_s == primary.dt_s.min()], "epsilon",
               rf"$\Delta t={scientific(primary.dt_s.min())}$ s")]
    for row, (subset, x, subtitle) in enumerate(slices):
        for column, metric in enumerate(["R16", "R26"]):
            lines(axes[row, column], subset, x, metric, signed=True)
            axes[row, column].set_title(subtitle, fontsize=11)
            axes[row, column].set_ylabel(r"$R_{16}$ [m]" if metric == "R16" else r"$R_{26}$")
            axes[row, column].set_xlabel(r"Time step $\Delta t$ [s]" if x == "dt_s"
                                          else r"Starting amplitude $\epsilon$")
    figure.suptitle("DBA: signed dispersion entries (analytic values are zero)", fontsize=14)
    figure.get_layout_engine().set(rect=(0, 0.045, 1, 0.845))
    legend(figure)
    figure.text(0.5, 0.014, r"Symmetric-log vertical axes, linear within $\pm10^{-9}$; dotted line = zero.",
                ha="center", fontsize=9)
    save(figure, output, "signed-dispersion")


def time_difference(primary: pd.DataFrame, output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(9.3, 4.5), layout="constrained")
    metric = "difference_from_finest_dt"
    for axis, epsilon in zip(axes, [1e-2, 1e-3]):
        subset = primary[primary.epsilon == epsilon]
        lines(axis, subset, "dt_s", metric)
        axis.set_xlabel(r"Time step $\Delta t$ [s]")
        axis.set_ylabel(r"$\max_{ij}|A_{ij}(\Delta t)-A_{ij}(10^{-13}\,\mathrm{s})|$")
        axis.set_title(rf"$\epsilon={scientific(epsilon)}$", fontsize=11)
        coarse = subset[(subset.richardson_levels == 1) & (subset.dt_s >= 3e-12)].sort_values("dt_s")
        if len(coarse) >= 2:
            x = coarse.dt_s.to_numpy()
            y = coarse[metric].iloc[-1] * (x / x[-1])**2
            axis.plot(x, y, ":", color="0.3", lw=1.4)
            axis.annotate(r"$\propto\Delta t^2$", (x[0], y[0]), xytext=(-5, 13),
                          textcoords="offset points", fontsize=10)
    figure.suptitle("DBA: time-step differences at fixed differentiation settings", fontsize=14)
    figure.get_layout_engine().set(rect=(0, 0.10, 1, 0.75))
    legend(figure)
    figure.text(0.5, 0.025, "The finest numerical map is not an exact reference; "
                "its zero self-difference is omitted.", ha="center", fontsize=9)
    save(figure, output, "time-differences")


def make_plots(frame: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "svg.fonttype": "none", "axes.spines.top": False,
                         "axes.spines.right": False})
    primary = frame[(frame.ranks == 1) & (frame.status == "OK")]
    cuts(primary, output, "dt_s", "epsilon", [1e-2, 1e-3], "errors-vs-dt")
    cuts(primary, output, "epsilon", "dt_s", [1e-11, 1e-13], "errors-vs-epsilon")
    heatmaps(frame, output)
    best_vs_dt(primary, output)
    signed_entries(primary, output)
    time_difference(primary, output)
