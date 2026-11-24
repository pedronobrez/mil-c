"""
What `qc.check` found, on a page: the flow rate over time and, under it,
the median of a few populations in each bin of events, its slow level and
the short stretches that stood apart from their neighbours, with the
stretches to remove shaded. The same series the window's dialog draws
(`ui/quality_dialog.py`), for a file instead of a screen.

No Qt here, and matplotlib is imported only when a page is drawn.
"""

from __future__ import annotations

import os

import numpy as np

from . import qc

#: the page, in inches: the flow rate, and each population under it
WIDTH, RATE_HEIGHT, SIGNAL_HEIGHT = 8.0, 2.2, 1.8
SHADE, LINE, SLOW, LIMIT = "#d1495b", "#2f6690", "#1b2a38", "#888888"
FORMATS = ("pdf", "png")


def timestep(keywords: dict) -> float:
    """Seconds in one unit of Time, as FCS 3 says it; 0 where the file
    does not."""
    try:
        return float(keywords.get("$TIMESTEP", "") or 0)
    except ValueError:
        return 0.0


def stretches(intervals, step: float = 0.0, most: int = 4) -> str:
    """Time intervals as text: in seconds where the file says how long a
    unit of Time is, otherwise in its units."""
    scale, unit = (step, " s") if step else (1.0, "")
    said = [f"{lo * scale:.4g}–{hi * scale:.4g}{unit}" for lo, hi in intervals[:most]]
    if len(intervals) > most:
        said.append(f"and {len(intervals) - most} more")
    return ", ".join(said)


def worst(result: qc.Result, most: int = 3) -> list[int]:
    """The populations to draw: those suspect for longest first, then the
    others in the file's order."""
    b = result.bins
    bad = (b.bad_signal | b.bad_drift | b.bad_shift).sum(axis=0)
    return [int(j) for j in np.argsort(-bad, kind="stable")[:most]]


def series(result: qc.Result, j: int, step: float = 0.0) -> dict:
    """What the plots of population `j` are drawn from, in the units of
    their axes: the time of each bin, its flow rate, the population's
    median (not a number where a bin has none of it), the limits of both,
    the slow level the drift is judged on with its own limits, the short
    stretches that stood apart — (first time, last time, their level,
    their neighbours' level) — and the stretches shaded."""
    b = result.bins
    step = step or 1.0
    scale = result.scales[b.channel[j]]
    back = scale.inverse if scale is not None else np.asarray
    median = np.where(np.isfinite(b.medians[:, j]), back(b.medians[:, j]), np.nan)
    slow = np.where(np.isfinite(b.slow[:, j]), back(b.slow[:, j]), np.nan)
    return {"time": (b.start + b.end) / 2 * step, "rate": b.rate / step,
            "rate_limits": [v / step for v in b.rate_limits],
            "median": median, "limits": np.asarray(back(b.limits[j])).tolist(),
            "bad_rate": b.bad_rate, "bad_signal": b.bad_signal[:, j],
            "slow": slow, "slow_limits": np.asarray(back(b.slow_limits[j])).tolist(),
            "bad_drift": b.bad_drift[:, j], "bad_shift": b.bad_shift[:, j],
            "shifts": [(float(b.start[a]) * step, float(b.end[e - 1]) * step,
                        float(back(level)), float(back(near)))
                       for k, a, e, level, near in b.shifts if k == j],
            "shaded": [(lo * step, hi * step) for lo, hi in result.intervals]}


def figure(name: str, result: qc.Result, step: float = 0.0, signals: int = 3):
    """One sample's page, as a matplotlib figure."""
    from matplotlib.figure import Figure
    if not result.checked:
        fig = Figure(figsize=(WIDTH, RATE_HEIGHT))
        fig.text(0.06, 0.62, name, fontsize=12, fontweight="bold")
        fig.text(0.06, 0.38, result.summary, fontsize=9, wrap=True)
        return fig
    chosen = worst(result, signals)
    heights = [RATE_HEIGHT] + [SIGNAL_HEIGHT] * len(chosen)
    fig = Figure(figsize=(WIDTH, sum(heights) + 1.0), layout="constrained")
    axes = np.atleast_1d(fig.subplots(len(heights), 1, sharex=True, height_ratios=heights))
    fig.suptitle(f"{name}\n{result.summary}", fontsize=9, x=0.02, ha="left", wrap=True)
    rows = [(axes[0], "rate", "bad_rate", "rate_limits",
             "Events per second" if step else "Events per unit of time")]
    rows += [(ax, "median", "bad_signal", "limits", result.populations[j])
             for ax, j in zip(axes[1:], chosen)]
    for k, (ax, y, bad, limits, label) in enumerate(rows):
        s = series(result, chosen[k - 1] if k else 0, step)
        for lo, hi in s["shaded"]:
            ax.axvspan(lo, hi, color=SHADE, alpha=0.18, linewidth=0)
        for level in s[limits]:
            ax.axhline(level, color=LIMIT, linewidth=0.8, linestyle="--")
        ax.plot(s["time"], s[y], color=LINE, linewidth=0.9)
        marked = s[bad] & np.isfinite(s[y])
        ax.plot(s["time"][marked], s[y][marked], "o", color=SHADE, markersize=2.5)
        if k:
            for level in s["slow_limits"]:
                ax.axhline(level, color=LIMIT, linewidth=0.8, linestyle=":")
            ax.plot(s["time"], s["slow"], color=SLOW, linewidth=1.4)
            ax.plot(s["time"], np.where(s["bad_drift"], s["slow"], np.nan), color=SHADE,
                    linewidth=2.2)
            for lo, hi, level, near in s["shifts"]:
                ax.plot([lo, hi], [near, near], color=SLOW, linewidth=1.0, linestyle="--")
                ax.plot([lo, hi], [level, level], color=SHADE, linewidth=2.6)
        if not k:
            ax.set_yscale("log")   # its limits are a factor either side of the usual rate
        ax.set_ylabel(label, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.margins(x=0.01)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[-1].set_xlabel("Time (s)" if step else "Time (as the file stores it)", fontsize=8)
    return fig


def save(path: str, pages: list[tuple[str, qc.Result, float]], fmt: str = "png",
         signals: int = 3, dpi: int = 150) -> list[str]:
    """The pages of (name, result, timestep): one PDF of them all where
    `path` ends in .pdf, otherwise a folder with a file for each sample.
    The files written."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    if path.lower().endswith(".pdf"):
        from matplotlib.backends.backend_pdf import PdfPages
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with PdfPages(path) as pdf:
            for name, result, step in pages:
                pdf.savefig(figure(name, result, step, signals))
        return [path]
    if fmt not in FORMATS:
        raise ValueError(f"cannot write .{fmt}: use " + ", ".join("." + f for f in FORMATS))
    os.makedirs(path, exist_ok=True)
    written = []
    for name, result, step in pages:
        fig = figure(name, result, step, signals)
        FigureCanvasAgg(fig)
        written.append(os.path.join(path, f"{name}_qc.{fmt}"))
        fig.savefig(written[-1], format=fmt, dpi=dpi)
    return written
