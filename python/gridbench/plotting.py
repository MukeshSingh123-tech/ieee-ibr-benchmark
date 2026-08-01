"""Figure generation for the report.

House style, applied once here so every figure in the report reads as one set:

  * ONE axis per figure. Never a dual y-axis -- two measures of different scale
    become two panels or an indexed common base.
  * Categorical colours assigned in FIXED slot order and never cycled. A series
    keeps its colour when other series are filtered out, so colour tracks the
    entity rather than its rank.
  * Thin marks, recessive grid and spines, legend whenever there are two or more
    series (a single series is named by the title instead).
  * Text uses ink colours, never the series colour. A coloured marker beside a
    label carries the identity.

Output is print-oriented PNG at 200 dpi for the LaTeX report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402
import numpy as np                     # noqa: E402

from .config import FIGURE_DIR         # noqa: E402

# --- validated categorical palette, fixed order -----------------------------
SERIES = [
    "#2a78d6",   # 1 blue
    "#eb6834",   # 2 orange
    "#1baf7a",   # 3 aqua
    "#eda100",   # 4 yellow
    "#e87ba4",   # 5 magenta
    "#008300",   # 6 green
    "#4a3aa7",   # 7 violet
    "#e34948",   # 8 red
]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e3e2df"

# status colours are RESERVED -- never reused as a categorical series
STATUS_GOOD = "#008300"
STATUS_BAD = "#e34948"


def apply_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK,
        "axes.titlesize": 11,
        "axes.titleweight": "600",
        "axes.titlecolor": INK,
        "axes.labelsize": 9.5,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.7,
        "xtick.color": INK_SOFT,
        "ytick.color": INK_SOFT,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.labelcolor": INK_SOFT,
        "lines.linewidth": 2.0,
        "lines.markersize": 5.0,
        "font.family": "DejaVu Sans",
        "figure.dpi": 120,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })


def _finish(ax, title: str, xlabel: str, ylabel: str, note: str | None = None) -> None:
    ax.set_title(title, loc="left", pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    if note:
        ax.text(0.0, -0.20, note, transform=ax.transAxes, fontsize=7.5,
                color=INK_SOFT, va="top", ha="left", wrap=True)


def save(fig, name: str) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURE_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# =============================================================================
# figure builders
# =============================================================================
def line_by_series(
    x: Sequence[float],
    series: dict[str, Sequence[float]],
    title: str,
    xlabel: str,
    ylabel: str,
    note: str | None = None,
    threshold: tuple[float, str] | None = None,
    figsize: tuple[float, float] = (6.4, 4.0),
    logy: bool = False,
    annotate_last: bool = True,
):
    """Line chart, one line per series, colours in fixed slot order.

    `threshold` draws a reference line (e.g. the SCR >= 3.0 interconnection
    screen) in ink rather than a series colour, so it reads as an annotation and
    not as another data series.
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    for i, (label, ys) in enumerate(series.items()):
        colour = SERIES[i % len(SERIES)]
        ys = np.asarray(ys, dtype=float)
        ax.plot(x, ys, color=colour, marker="o", label=label,
                markeredgecolor=SURFACE, markeredgewidth=1.0)
        # selective direct label on the last finite point only
        if annotate_last and len(series) <= 4:
            finite = np.flatnonzero(np.isfinite(ys))
            if finite.size:
                k = finite[-1]
                ax.annotate(f"{ys[k]:.4g}", (x[k], ys[k]), textcoords="offset points",
                            xytext=(6, 0), fontsize=8, color=INK_SOFT, va="center")

    if threshold is not None:
        value, text = threshold
        ax.axhline(value, color=INK_SOFT, linewidth=1.0, linestyle=(0, (4, 3)))
        ax.text(x[0], value, f"  {text}", fontsize=7.5, color=INK_SOFT,
                va="bottom", ha="left")

    if logy:
        ax.set_yscale("log")
    if len(series) >= 2:
        ax.legend(loc="best")
    _finish(ax, title, xlabel, ylabel, note)
    return fig


def grouped_bars(
    categories: Sequence[str],
    series: dict[str, Sequence[float]],
    title: str,
    xlabel: str,
    ylabel: str,
    note: str | None = None,
    logy: bool = False,
    figsize: tuple[float, float] = (6.4, 4.0),
    value_fmt: str = "{:.0f}",
):
    """Grouped bars with a 2px surface gap between adjacent fills."""
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    n_groups = len(categories)
    n_series = len(series)
    slot = 0.8 / n_series
    idx = np.arange(n_groups)

    for i, (label, ys) in enumerate(series.items()):
        colour = SERIES[i % len(SERIES)]
        pos = idx - 0.4 + slot * (i + 0.5)
        ys = np.asarray(ys, dtype=float)
        ax.bar(pos, ys, width=slot * 0.88, color=colour, label=label,
               edgecolor=SURFACE, linewidth=1.5)
        for p, v in zip(pos, ys):
            if np.isfinite(v):
                ax.annotate(value_fmt.format(v), (p, v), textcoords="offset points",
                            xytext=(0, 3), ha="center", fontsize=7.2, color=INK_SOFT)

    ax.set_xticks(idx)
    ax.set_xticklabels(categories)
    if logy:
        ax.set_yscale("log")
    else:
        # headroom so the value labels above the tallest bar are not clipped
        top = np.nanmax([np.nanmax(np.asarray(v, dtype=float)) for v in series.values()])
        if np.isfinite(top) and top > 0:
            ax.set_ylim(0, top * 1.15)

    _finish(ax, title, xlabel, ylabel, note)
    if n_series >= 2:
        # Above the axes, in one row. A frameless legend placed "best" inside a
        # grouped-bar chart sits on top of the bars, and the bars show through
        # the swatches -- which misreports the series colours.
        # No mode="expand": it forces equal-width columns and truncates labels.
        ax.legend(loc="upper left", bbox_to_anchor=(0, 1.13),
                  ncols=min(n_series, 4), borderaxespad=0,
                  columnspacing=1.4, handletextpad=0.5)
        ax.set_title(title, loc="left", pad=32)
    return fig


def heatmap(
    values: np.ndarray,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    title: str,
    xlabel: str,
    ylabel: str,
    cbar_label: str,
    note: str | None = None,
    figsize: tuple[float, float] = (6.6, 4.2),
    fmt: str = "{:.3f}",
):
    """Sequential heatmap: ONE hue, light -> dark. Never a rainbow.

    Magnitude is the job here, so a single-hue ramp is the correct encoding --
    a rainbow would imply category boundaries that do not exist in the data.
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "blues", ["#f2f6fc", "#2a78d6", "#12305c"])
    masked = np.ma.masked_invalid(np.asarray(values, dtype=float))
    cmap.set_bad("#eeeeec")

    im = ax.imshow(masked, cmap=cmap, aspect="auto")
    ax.set_xticks(np.arange(len(col_labels)), labels=col_labels)
    ax.set_yticks(np.arange(len(row_labels)), labels=row_labels)

    # 2px surface gridlines between cells
    ax.set_xticks(np.arange(len(col_labels) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(row_labels) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2.0)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="minor", length=0)

    vmax = np.nanmax(masked) if masked.count() else 1.0
    for i in range(masked.shape[0]):
        for j in range(masked.shape[1]):
            v = masked[i, j]
            if v is np.ma.masked:
                ax.text(j, i, "n/a", ha="center", va="center",
                        fontsize=7, color=INK_SOFT)
            else:
                ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=7.2,
                        color="#ffffff" if v > 0.6 * vmax else INK)

    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label(cbar_label, fontsize=8.5, color=INK_SOFT)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, color=GRID)

    ax.set_title(title, loc="left", pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for s in ax.spines.values():
        s.set_visible(False)
    if note:
        ax.text(0.0, -0.22, note, transform=ax.transAxes, fontsize=7.5,
                color=INK_SOFT, va="top", ha="left")
    return fig
