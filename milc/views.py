"""
Three more ways of looking at a gated population, each answering a
question the plot of two parameters does not:

- violins and boxes: how one marker is distributed in each sample, or in
  each group of them, side by side;
- the matrix of every pair of parameters: where spillover and an
  unexpected correlation show, without choosing the pair beforehand;
- backgating: where the events of a population fall on each gate above
  it — whether the scatter gate cut into the cells one ended up with.

Like `density`, everything here works in the plot's own coordinates,
after the display scale. The numbers of a box are therefore taken on the
axis it is drawn on: quartiles are the same on any monotonic scale (to
within the interpolation between two events), the reach of the whiskers
is not, and it is the drawn one that is meant.

The contents and their matplotlib drawing live here rather than in
`figures`, and are registered with it when this module is imported.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import density, figures, order
from .gating import (BooleanGate, EllipseGate, GateNode, PolygonGate, QuadrantGate,
                     RectangleGate)
from .transforms import decade_ticks, same

VIOLIN, BOX = "violin", "box"

#: Tukey's rule: a whisker ends at the furthest event within this many
#: interquartile ranges of the box; what lies beyond is an outlier
WHISKER = 1.5
#: outliers drawn for one box. One sample of 200,000 events has thousands
#: beyond the whiskers, and a figure with every one is a black bar
MAX_OUTLIERS = 200

#: bins along each parameter of the matrix. A cell of a 10 x 10 grid is a
#: few centimetres on any screen or page: more bins than this are not seen
MATRIX_BINS = 64
#: 12 parameters are 144 cells; beyond that a cell is too small to read
MAX_PARAMETERS = 12

BACKGATE_BINS = 96
#: events of the population drawn as dots on each step of a backgate
MAX_HIGHLIGHTED = 5000

#: half the room of its row a violin's widest point takes, and a box
_HALF, _BOX_HALF = 0.42, 0.26


# --- violins and boxes -----------------------------------------------------------

@dataclass
class Summary:
    """The numbers of a box plot. `low` and `high` are where the whiskers
    end: at events, not at the limits of the rule."""
    n: int
    median: float
    q1: float
    q3: float
    low: float
    high: float
    outliers: int


def summary(x: np.ndarray) -> Summary:
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0:
        return Summary(0, *[float("nan")] * 5, 0)
    q1, median, q3 = np.percentile(x, [25, 50, 75])
    reach = WHISKER * (q3 - q1)
    within = x[(x >= q1 - reach) & (x <= q3 + reach)]
    return Summary(int(x.size), float(median), float(q1), float(q3),
                   float(within.min()), float(within.max()), int(x.size - within.size))


@dataclass
class Distribution:
    """One sample, or one group of them pooled, as a violin or a box."""
    label: str
    #: what its colour stands for, as on an overlay
    group: str
    summary: Summary
    #: the violin: `density.histogram`'s own centres and heights, so a
    #: violin and the overlaid histogram of the same events are one curve
    x: np.ndarray
    density: np.ndarray
    #: the outliers that are drawn: all, or `MAX_OUTLIERS` spread evenly
    #: through them in order, the two extremes always among them
    outliers: np.ndarray


def distribution(x: np.ndarray, lo: float, hi: float, label: str, group: str = "",
                 max_outliers: int = MAX_OUTLIERS) -> Distribution:
    """`x` on the display scale, [lo, hi] the ends of its axis."""
    x = np.asarray(x, dtype=np.float64)
    s = summary(x)
    centres, heights = density.histogram(x, lo, hi)
    far = np.sort(x[(x < s.low) | (x > s.high)])
    if far.size > max_outliers:
        far = far[np.round(np.linspace(0, far.size - 1, max_outliers)).astype(int)]
    return Distribution(label, group or label, s, centres, heights, far)


def violin(d: Distribution) -> tuple[np.ndarray, np.ndarray]:
    """(positions, half-widths in rows) of a violin's outline, between the
    first and the last place it is wider than a hair: the tails of a
    smoothed histogram reach the ends of the axis as a line of no width."""
    wide = np.flatnonzero(d.density > 0.5)
    if wide.size == 0:
        return d.x[:0], d.density[:0]
    keep = slice(wide[0], wide[-1] + 1)
    return d.x[keep], _HALF * d.density[keep] / 100.0


@dataclass
class DistributionContent:
    series: list[Distribution]
    kind: str = VIOLIN
    xlabel: str = ""
    xticks: list | None = None
    xrange: tuple | None = None

    def groups(self) -> list[str]:
        return order.ordered(s.group for s in self.series)

    def category_names(self) -> list[tuple[str, str]]:
        names = self.groups()
        return list(zip(names, figures.categorical(len(names))))

    def note(self) -> str:
        """What the marks stand for: said wherever the plot is shown."""
        if self.kind == VIOLIN:
            return "dot: median; bar: quartiles"
        total = sum(s.summary.outliers for s in self.series)
        drawn = sum(s.outliers.size for s in self.series)
        capped = (f"; {drawn:,d} of {total:,d} outliers drawn, at most {MAX_OUTLIERS} a box"
                  if drawn < total else "")
        return (f"box: quartiles and median; whiskers: the furthest event within "
                f"{WHISKER:g} IQR of the box{capped}")


def _draw_distributions(fig, c: DistributionContent, style) -> None:
    from matplotlib.colors import to_rgba
    from matplotlib.patches import Patch, Rectangle
    ax = fig.subplots()
    names = c.category_names()
    colour = {g: style.color(g, default) for g, default in names}
    for k, s in enumerate(c.series):
        m, col = s.summary, colour[s.group]
        if not m.n:
            continue
        if c.kind == VIOLIN:
            x, half = violin(s)
            ax.fill_between(x, k - half, k + half, facecolor=to_rgba(col, 0.35),
                            edgecolor=col, linewidth=0.9)
            ax.plot([m.q1, m.q3], [k, k], color="black", linewidth=2.2, solid_capstyle="butt")
            ax.plot([m.median], [k], "o", markersize=3.4, markerfacecolor="white",
                    markeredgecolor="black", markeredgewidth=0.6)
            continue
        ax.add_patch(Rectangle((m.q1, k - _BOX_HALF), m.q3 - m.q1, 2 * _BOX_HALF,
                               facecolor=to_rgba(col, 0.35), edgecolor=col, linewidth=0.9))
        ax.plot([m.median, m.median], [k - _BOX_HALF, k + _BOX_HALF], color="black",
                linewidth=1.1, solid_capstyle="butt")
        cap = _BOX_HALF / 2
        for a, b in ((m.low, m.q1), (m.q3, m.high)):
            ax.plot([a, b], [k, k], color=col, linewidth=0.9)
        for end in (m.low, m.high):
            ax.plot([end, end], [k - cap, k + cap], color=col, linewidth=0.9)
        ax.scatter(s.outliers, np.full(s.outliers.size, float(k)),
                   s=(style.point_size * 0.6) ** 2, color=col, linewidths=0, alpha=0.7)
    ax.set_yticks(range(len(c.series)), [style.label(s.label) for s in c.series])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(c.series) - 0.4)
    if c.xticks:
        ax.set_xticks([p for p, _ in c.xticks], [figures._tick_text(t) for _, t in c.xticks])
    if c.xrange:
        ax.set_xlim(*c.xrange)
    ax.set_xlabel(c.xlabel if style.xlabel is None else style.xlabel)
    ax.set_ylabel("" if style.ylabel is None else style.ylabel)
    figures._frame(ax, style)
    # one colour per row: the axis already names each
    if len(names) != len(c.series):
        handles = [Patch(facecolor=to_rgba(colour[g], 0.35), edgecolor=colour[g],
                         label=style.label(g)) for g, _ in names]
        figures._legend(fig, [ax], handles, style, "")


# --- every pair of parameters -----------------------------------------------------

def _bins(x: np.ndarray, lo: float, hi: float, bins: int) -> np.ndarray:
    """The bin of each event; those beyond an end of the axis go in the
    bin at that end, as `density.histogram` has them."""
    at = (np.asarray(x, dtype=np.float64) - lo) * (bins / (hi - lo))
    return np.clip(at, 0, bins - 1).astype(np.intp)


def _count2(a: np.ndarray, b: np.ndarray, bins: int) -> np.ndarray:
    return np.bincount(a * bins + b, minlength=bins * bins).reshape(bins, bins)


def pairs(columns: list[np.ndarray], ranges: list[tuple[float, float]],
          bins: int = MATRIX_BINS) -> np.ndarray:
    """k parameters of the same events, binned two by two: k x k x bins x
    bins, where `[i, j, a, b]` counts the events in bin `a` of parameter i
    and bin `b` of parameter j. `[j, i]` is the transpose of `[i, j]`, and
    `[i, i]` holds the histogram of parameter i on its diagonal. Every
    cell sums to the number of events.

    Each event is put in its bin once per parameter and each pair is then
    one `bincount`; the result is 16 KiB a cell at 64 bins whatever the
    number of events."""
    k = len(columns)
    index = [_bins(x, lo, hi, bins) for x, (lo, hi) in zip(columns, ranges)]
    out = np.zeros((k, k, bins, bins), dtype=np.uint32)
    for i in range(k):
        for j in range(i, k):
            out[i, j] = _count2(index[i], index[j], bins)
            out[j, i] = out[i, j].T
    return out


def plain_ticks(lo: float, hi: float) -> list[tuple[float, str]]:
    """A few round numbers for an axis with no scale of its own, short
    enough for a cell of the matrix: `100K`."""
    from matplotlib.ticker import MaxNLocator
    at = [v for v in MaxNLocator(3).tick_values(lo, hi) if lo <= v <= hi]
    return [(float(v), f"{v / 1000:g}K" if abs(v) >= 1000 else f"{v:g}") for v in at]


@dataclass
class MatrixContent:
    """Cell (i, j), row i from the top, is parameter j across against
    parameter i up; on the diagonal, each parameter's histogram."""
    labels: list[str]
    counts: np.ndarray
    ranges: list[tuple[float, float]]
    #: per parameter, [(position, label)] on its axis
    ticks: list[list]
    uses_colormap = True

    def category_names(self) -> list[tuple[str, str]]:
        return []


def _draw_matrix(fig, c: MatrixContent, style) -> None:
    """One pair of axes for the whole grid, a cell a unit square. Every
    density cell is part of one bitmap (`imshow`), also in SVG and PDF:
    the file weighs the same for a hundred events or a million, and the
    text and the histograms stay vector."""
    k, bins = len(c.labels), c.counts.shape[-1]
    cmap = figures.colormap(style.colormap, style.colormap_colors)
    ax = fig.subplots()
    image = np.zeros((k * bins, k * bins, 4))
    centres = (np.arange(bins) + 0.5) / bins
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            # each cell on its own logarithmic ramp: a pair of dim markers
            # is as legible as a pair of bright ones
            cell = np.log1p(c.counts[i, j][::-1].astype(np.float64))
            rgba = cmap(cell / cell.max() if cell.max() > 0 else cell)
            rgba[cell == 0] = 0.0
            image[i * bins:(i + 1) * bins, j * bins:(j + 1) * bins] = rgba
    # "none": a vector file carries the bins themselves, not the bins
    # sampled again at the figure's resolution
    ax.imshow(image, extent=(0, k, 0, k), origin="upper", aspect="auto",
              interpolation="none")
    for i, label in enumerate(c.labels):
        base = k - 1 - i
        h = np.diagonal(c.counts[i, i]).astype(np.float64)
        top = h.max() if h.max() > 0 else 1.0
        # the upper quarter of the cell is the name's
        ax.fill_between(i + centres, base, base + 0.72 * h / top, step="mid",
                        color="#6b7280", linewidth=0)
        ax.text(i + 0.5, base + 0.96, style.label(label.split()[0] if label else label),
                ha="center", va="top", fontsize=style.font_size, clip_on=True)
    for at in range(1, k):
        ax.axhline(at, color="#c4c8d0", linewidth=0.4)
        ax.axvline(at, color="#c4c8d0", linewidth=0.4)
    where, text = [], []
    for i, ((lo, hi), ticks) in enumerate(zip(c.ranges, c.ticks)):
        for p, t in ticks:
            if lo <= p <= hi:
                where.append((i + (p - lo) / (hi - lo), k - 1 - i + (p - lo) / (hi - lo)))
                text.append(figures._tick_text(t))
    # the decades are marked on every column and row, and written only
    # where a cell is wide enough for them not to run into each other
    roomy = min(style.width_mm, style.height_mm) / k >= 20
    labels = text if roomy else [""] * len(text)
    ax.set_xticks([x for x, _ in where], labels, fontsize=style.font_size - 2)
    ax.set_yticks([y for _, y in where], labels, fontsize=style.font_size - 2)
    ax.set_xlim(0, k)
    ax.set_ylim(0, k)
    ax.set_xlabel("" if style.xlabel is None else style.xlabel)
    ax.set_ylabel("" if style.ylabel is None else style.ylabel)
    figures._frame(ax, style)


# --- backgating -------------------------------------------------------------------

@dataclass
class Step:
    """One gate on the way down to a population."""
    name: str
    path: str
    gate: object
    #: the parameters it is looked at on: the gate's own; for a
    #: combination, those of the first population it is made of that has
    #: any; none when there is nothing to draw it on
    dims: list[str]
    #: the events the gate was applied to: its parent's
    parent: np.ndarray


def _dims(node: GateNode) -> list[str]:
    if not isinstance(node.gate, BooleanGate):
        return list(node.gate.dims)
    for sister in node.parent.children:
        if sister.name in node.gate.operands and not isinstance(sister.gate, BooleanGate):
            return list(sister.gate.dims)
    return []


def steps(tree: GateNode, path: str, masks: dict[str, np.ndarray]) -> list[Step]:
    """Every gate from the top of the tree down to the population at
    `path`, its own included — the last step shows the gate that made it.
    `masks` is what `tree.apply` gave."""
    chain, node = [], tree.find(path)
    while node.gate is not None:
        chain.append(node)
        node = node.parent
    return [Step(n.name, n.path, n.gate, _dims(n), masks[n.parent.path])
            for n in reversed(chain)]


def _shown(values, scale) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return scale(values) if scale is not None else values


def _carried(points: np.ndarray, drawn: list, shown: list) -> np.ndarray:
    """Points in the coordinates a gate was drawn in, where they fall on
    axes that may have other scales by now."""
    out = np.array(points, dtype=np.float64)
    for k in range(2):
        if not same(drawn[k], shown[k]):
            stored = drawn[k].inverse(out[:, k]) if drawn[k] is not None else out[:, k]
            out[:, k] = _shown(stored, shown[k])
    return out


def outline(gate, scale_of, ranges: list[tuple[float, float]]) -> list[np.ndarray]:
    """A gate as lines (each n x 2) on axes drawn with `scale_of(dim)`,
    whose ends are `ranges`, one for each of the gate's dimensions.

    A rectangle's open side is drawn at the end of its axis. A gate of
    one parameter is its bounds as vertical lines from 0 to 1: of whatever
    height the plot has. A polygon or an ellipse is straight, or round,
    on the scales it is tested on; shown on others, it is the curve it
    has become. A quadrant is the two lines of its crossing. A
    combination has no outline."""
    if isinstance(gate, BooleanGate):
        return []
    shown = [scale_of(d) for d in gate.dims]
    if isinstance(gate, RectangleGate):
        ends = [(r[0] if lo is None else float(_shown(lo, s)),
                 r[1] if hi is None else float(_shown(hi, s)))
                for lo, hi, s, r in zip(gate.mins, gate.maxs, shown, ranges)]
        if len(gate.dims) == 1:
            return [np.array([[float(_shown(v, shown[0])), 0.0],
                              [float(_shown(v, shown[0])), 1.0]])
                    for v in (gate.mins[0], gate.maxs[0]) if v is not None]
        (x0, x1), (y0, y1) = ends
        return [np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]])]
    if isinstance(gate, QuadrantGate):
        px, py = float(_shown(gate.x, shown[0])), float(_shown(gate.y, shown[1]))
        return [np.array([[px, ranges[1][0]], [px, ranges[1][1]]]),
                np.array([[ranges[0][0], py], [ranges[0][1], py]])]
    drawn = list(gate.scales or [None, None])
    if isinstance(gate, PolygonGate):
        own = np.array(gate.vertices, dtype=np.float64)
        for k in range(2):
            own[:, k] = _shown(own[:, k], drawn[k])
        own = np.vstack([own, own[:1]])
        if not all(same(a, b) for a, b in zip(drawn, shown)):
            # a straight side is a curve on another scale: follow it
            t = np.linspace(0.0, 1.0, 24, endpoint=False)[:, None]
            own = np.vstack([a + t * (b - a) for a, b in zip(own[:-1], own[1:])] + [own[-1:]])
        return [_carried(own, drawn, shown)]
    if isinstance(gate, EllipseGate):
        t = np.linspace(0, 2 * np.pi, 121)
        c, s = np.cos(gate.angle), np.sin(gate.angle)
        u, v = gate.radii[0] * np.cos(t), gate.radii[1] * np.sin(t)
        own = np.column_stack([gate.center[0] + u * c - v * s, gate.center[1] + u * s + v * c])
        return [_carried(own, drawn, shown)]
    return []


@dataclass
class Panel:
    """One step, ready to draw."""
    title: str
    xlabel: str = ""
    ylabel: str = ""
    xticks: list | None = None
    yticks: list | None = None
    xrange: tuple | None = None
    yrange: tuple | None = None
    #: two parameters: the parent's events binned, `[y, x]`, and the
    #: population's own as points (n x 2, at most `MAX_HIGHLIGHTED`)
    counts: np.ndarray | None = None
    points: np.ndarray | None = None
    #: one parameter: centres, and the events per bin of the parent and of
    #: the population, smoothed alike so that one sits inside the other
    x: np.ndarray | None = None
    parent: np.ndarray | None = None
    chosen: np.ndarray | None = None
    outline: list = field(default_factory=list)


#: the first of the two things a backgate colours; the second is the
#: population by its own name
EVERYONE = "All of the parent's events"


@dataclass
class BackgateContent:
    population: str
    panels: list[Panel]
    #: events of the population, and how many of them each panel draws
    events: int = 0
    drawn: int = 0

    def category_names(self) -> list[tuple[str, str]]:
        return [(EVERYONE, "#6f7682"), (self.population, figures.CATEGORICAL[1])]


def backgate(tree: GateNode, path: str, masks: dict[str, np.ndarray], display, scale_of,
             range_of, label_of) -> BackgateContent:
    """The population at `path` on each gate above it. `display(dim)` is
    a column of every event on its display scale, `scale_of(dim)` that
    scale, `range_of(dim)` the ends of its axis, `label_of(dim)` its name."""
    chosen = masks[path]
    rows = np.flatnonzero(chosen)
    if rows.size > MAX_HIGHLIGHTED:
        rows = np.sort(np.random.default_rng(0).choice(rows, MAX_HIGHLIGHTED, replace=False))
    panels = []
    for step in steps(tree, path, masks):
        if not step.dims:
            panels.append(Panel(step.name))
            continue
        x, xr, xs = display(step.dims[0]), range_of(step.dims[0]), scale_of(step.dims[0])
        common = dict(xlabel=label_of(step.dims[0]), xrange=xr,
                      xticks=decade_ticks(xs) if xs is not None else None)
        if len(step.dims) == 1:
            smooth = density.auto_smooth(int(step.parent.sum()))
            cx, everyone = density.histogram(x[step.parent], *xr, smooth=smooth,
                                             normalise=density.COUNT)
            _, own = density.histogram(x[chosen], *xr, smooth=smooth,
                                       normalise=density.COUNT)
            panels.append(Panel(step.name, ylabel="Events", x=cx, parent=everyone, chosen=own,
                                outline=outline(step.gate, scale_of, [xr]), **common))
            continue
        y, yr, ys = display(step.dims[1]), range_of(step.dims[1]), scale_of(step.dims[1])
        counts = _count2(_bins(y[step.parent], *yr, BACKGATE_BINS),
                         _bins(x[step.parent], *xr, BACKGATE_BINS), BACKGATE_BINS)
        panels.append(Panel(step.name, ylabel=label_of(step.dims[1]), yrange=yr,
                            yticks=decade_ticks(ys) if ys is not None else None,
                            counts=counts, points=np.column_stack([x[rows], y[rows]]),
                            outline=outline(step.gate, scale_of, [xr, yr]), **common))
    return BackgateContent(tree.find(path).name, panels, int(chosen.sum()), int(rows.size))


def _draw_backgate(fig, c: BackgateContent, style) -> None:
    """The parent's events are a bitmap and so are the population's dots
    (`rasterized`), in SVG and PDF too: up to 5,000 circles a panel would
    make a file nobody opens. Outlines, axes and text stay vector."""
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    from matplotlib.patches import Patch
    (_, grey), (_, red) = c.category_names()
    grey, red = style.color(EVERYONE, grey), style.color(c.population, red)
    # a bin with one event in it is still seen: the ramp starts a quarter
    # of the way from white
    faint = tuple(1.0 - 0.25 * (1.0 - v) for v in to_rgb(grey))
    ramp = LinearSegmentedColormap.from_list("parent", [faint, grey])
    cols = min(len(c.panels), 4)
    rows = -(-len(c.panels) // cols)
    grid = fig.subplots(rows, cols, squeeze=False)
    for extra in grid.ravel()[len(c.panels):]:
        extra.set_visible(False)
    for ax, p in zip(grid.ravel(), c.panels):
        ax.set_title(style.label(p.title))
        if p.counts is not None:
            shade = np.log1p(p.counts.astype(np.float64))
            ax.imshow(np.ma.masked_equal(shade, 0.0), extent=(*p.xrange, *p.yrange),
                      origin="lower", aspect="auto", interpolation="none", cmap=ramp,
                      vmin=0.0, vmax=max(float(shade.max()), 1e-9))
            ax.scatter(p.points[:, 0], p.points[:, 1], s=(style.point_size * 0.8) ** 2,
                       color=red, linewidths=0, alpha=0.8, rasterized=True)
            for line in p.outline:
                ax.plot(line[:, 0], line[:, 1], color="black", linewidth=0.8)
            ax.set_ylim(*p.yrange)
        elif p.x is not None:
            top = max(float(p.parent.max()), 1.0) * 1.05
            ax.fill_between(p.x, 0, p.parent, color=grey, alpha=0.7, linewidth=0)
            ax.fill_between(p.x, 0, p.chosen, color=red, alpha=0.8, linewidth=0)
            for line in p.outline:
                ax.plot(line[:, 0], line[:, 1] * top, color="black", linewidth=0.8)
            ax.set_ylim(0, top)
        else:
            ax.text(0.5, 0.5, "a combination:\nno plot of its own", ha="center", va="center",
                    transform=ax.transAxes, fontsize=style.font_size - 1)
            ax.set_xticks([])
            ax.set_yticks([])
        if p.xticks:
            ax.set_xticks([q for q, _ in p.xticks], [figures._tick_text(t) for _, t in p.xticks])
        if p.yticks:
            ax.set_yticks([q for q, _ in p.yticks], [figures._tick_text(t) for _, t in p.yticks])
        if p.xrange:
            ax.set_xlim(*p.xrange)
        ax.set_xlabel(p.xlabel if style.xlabel is None else style.xlabel)
        ax.set_ylabel(p.ylabel if style.ylabel is None else style.ylabel)
        figures._frame(ax, style)
    count = f" ({c.events:,d})" if style.legend_counts else ""
    handles = [Patch(facecolor=grey, alpha=0.7, linewidth=0, label=style.label(EVERYONE)),
               figures._dot(red, style.label(c.population) + count)]
    figures._legend(fig, list(grid.ravel()[:len(c.panels)]), handles, style, "")


figures._DRAW.update({DistributionContent: _draw_distributions, MatrixContent: _draw_matrix,
                      BackgateContent: _draw_backgate})
