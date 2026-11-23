"""
Figures for paper: what a plot shows (a *content*), how it should look (a
`FigureStyle`), and one function that draws the two into a matplotlib
figure and writes it — SVG, PDF, PNG or TIFF, at the size in millimetres a
journal asks for.

The window's own plots are for looking; these are for printing. Both are
drawn from the same content objects, so what is exported is what was on
the screen, restyled. Text stays text in SVG and PDF, so a label can still
be edited in Illustrator or Inkscape.

No Qt here, and matplotlib is imported only when a figure is drawn.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields

import numpy as np

from . import cluster_stats

#: Paul Tol's qualitative colours: distinguishable with the common
#: colour-vision deficiencies
CATEGORICAL = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377",
               "#EE7733", "#009988", "#882255", "#999933", "#332288", "#BBBBBB"]

COLORMAPS = ["viridis", "magma", "inferno", "plasma", "cividis", "turbo",
             "Blues", "Purples", "Greens", "Reds", "Greys", "YlOrRd", "YlGnBu",
             "coolwarm", "RdBu_r", "PiYG"]

#: a ramp made of the user's own colours, low to high
CUSTOM = "custom"
DEFAULT_CUSTOM = ["#f7f7fb", "#8f7fd6", "#2a1a6e"]


def colormap(name: str = "viridis", colors: list[str] | None = None):
    """A matplotlib colour map: one of `COLORMAPS`, or `CUSTOM` running
    through `colors` (two or more, low to high). An unknown name is
    viridis, not an error: a preset may outlive the list."""
    from matplotlib import colormaps
    from matplotlib.colors import LinearSegmentedColormap
    if name == CUSTOM:
        colors = [c for c in (colors or []) if c] or DEFAULT_CUSTOM
        if len(colors) == 1:
            colors = ["#ffffff", colors[0]]
        return LinearSegmentedColormap.from_list("custom", colors)
    return colormaps[name] if name in colormaps else colormaps["viridis"]


def lookup_table(name: str = "viridis", colors: list[str] | None = None,
                 steps: int = 256) -> np.ndarray:
    """steps x 4, 0-255: the ramp for anything that is not matplotlib."""
    return (colormap(name, colors)(np.linspace(0.0, 1.0, steps)) * 255).astype(np.uint8)

FORMATS = {"svg": "SVG, vector (*.svg)", "pdf": "PDF, vector (*.pdf)",
           "png": "PNG (*.png)", "tiff": "TIFF (*.tiff)"}


def categorical(n: int) -> list[str]:
    return [CATEGORICAL[i % len(CATEGORICAL)] for i in range(n)]


# --- what a plot shows ---------------------------------------------------------

@dataclass
class Categories:
    """Each event in one of `names`, by index."""
    codes: np.ndarray
    names: list[str]
    colors: list[str]
    title: str = ""

    def counts(self) -> np.ndarray:
        return np.bincount(self.codes, minlength=len(self.names))


@dataclass
class Values:
    """Each event a number, coloured along a ramp between `lo` and `hi`."""
    values: np.ndarray
    title: str = ""
    lo: float | None = None
    hi: float | None = None
    #: [(value, label)] for the colour bar, where the values are on a
    #: scale of their own
    ticks: list | None = None

    def limits(self) -> tuple[float, float]:
        v = self.values[np.isfinite(self.values)]
        if v.size == 0:
            return 0.0, 1.0
        # the extreme percent on either side would spend the whole ramp
        lo = np.percentile(v, 1) if self.lo is None else self.lo
        hi = np.percentile(v, 99) if self.hi is None else self.hi
        return (float(lo), float(hi)) if hi > lo else (float(lo), float(lo) + 1.0)


def categories_of(labels, title: str = "", order: list[str] | None = None) -> Categories:
    """From per-event labels of any kind; numbers sort as numbers."""
    labels = np.asarray(labels)
    names, codes = np.unique(labels, return_inverse=True)
    names = [str(n) for n in names]
    if order is None:
        try:
            rank = np.argsort([float(n) for n in names], kind="stable")
        except ValueError:
            rank = np.arange(len(names))
    else:
        rank = np.array([names.index(n) for n in order if n in names]
                        + [i for i, n in enumerate(names) if n not in order])
    remap = np.empty(len(names), dtype=np.int64)
    remap[rank] = np.arange(len(names))
    names = [names[i] for i in rank]
    return Categories(remap[codes], names, categorical(len(names)), title)


def combine(a: Categories | None, b: Categories | None) -> Categories | None:
    """Two subdivisions as one: every combination that occurs, `Case · 3`,
    in the order of the first and then of the second."""
    if a is None or b is None:
        return a if b is None else b
    key = a.codes.astype(np.int64) * len(b.names) + b.codes
    present, codes = np.unique(key, return_inverse=True)
    names = [f"{a.names[k // len(b.names)] or '—'} · {b.names[k % len(b.names)] or '—'}"
             for k in present]
    return Categories(codes, names, categorical(len(names)),
                      " · ".join(t for t in (a.title, b.title) if t))


@dataclass
class MapContent:
    x: np.ndarray
    y: np.ndarray
    xlabel: str = ""
    ylabel: str = ""
    color: Categories | Values | None = None
    split: Categories | None = None
    #: [(position, label)] where an axis is not a plain number: a marker on
    #: a biexponential scale, ticked in decades
    xticks: list | None = None
    yticks: list | None = None
    #: panels per row when split; None lets the number of panels decide.
    #: Split by two columns, a row is one value of the first.
    columns: int | None = None


@dataclass
class HeatmapContent:
    table: cluster_stats.Table
    per_column: bool = True


@dataclass
class FrequencyContent:
    freq: cluster_stats.Table
    #: the group of each column of `freq`: what is compared, by colour
    groups: list[str]
    #: optionally, a second subdivision of the samples: one panel for each
    #: of its values — groups compared within each day, say
    facets: list[str] | None = None


@dataclass
class TreeContent:
    tree: dict
    #: one of `tree["markers"]`, or None for the metacluster
    color_by: str | None = None


@dataclass
class Series:
    """One sample (or one group of them) on an overlay."""
    label: str
    #: what its colour stands for: its own name, or the value of a column
    group: str
    #: a histogram: centres and heights
    x: np.ndarray | None = None
    y: np.ndarray | None = None
    #: contours: [(share of events inside, [n x 2 lines])], outermost first
    levels: list | None = None
    #: events outside the outermost contour, n x 2
    outliers: np.ndarray | None = None


@dataclass
class OverlayContent:
    """Several samples on one plot: histograms of one parameter, or the
    contours of two."""
    series: list[Series]
    xlabel: str = ""
    ylabel: str = ""
    xticks: list | None = None
    yticks: list | None = None
    xrange: tuple | None = None
    yrange: tuple | None = None
    #: histograms one above the other instead of on one baseline
    offset: bool = False

    @property
    def is_histogram(self) -> bool:
        return bool(self.series) and self.series[0].x is not None

    def groups(self) -> list[str]:
        return list(dict.fromkeys(s.group for s in self.series))


#: how far up each staggered histogram starts, as a share of the tallest
STAGGER = 0.62


def category_names(content) -> list[tuple[str, str]]:
    """(name, default colour) of what the content colours by category —
    what a person may want to rename or recolour."""
    if isinstance(content, OverlayContent):
        names = content.groups()
        return list(zip(names, categorical(len(names))))
    if isinstance(content, MapContent) and isinstance(content.color, Categories):
        return list(zip(content.color.names, content.color.colors))
    if isinstance(content, FrequencyContent):
        names = cluster_stats.ordered(content.groups)
        return list(zip(names, categorical(len(names))))
    if isinstance(content, TreeContent) and content.color_by is None:
        t = content.tree
        ids = sorted({int(c) for c in t["cluster"][t["sizes"] > 0]})
        return list(zip([str(i) for i in ids], categorical(len(ids))))
    return []


def uses_colormap(content) -> bool:
    return (isinstance(content, HeatmapContent)
            or (isinstance(content, MapContent) and isinstance(content.color, Values))
            or (isinstance(content, TreeContent) and content.color_by is not None))


# --- how it should look ----------------------------------------------------------

@dataclass
class FigureStyle:
    width_mm: float = 120.0
    height_mm: float = 90.0
    #: for PNG and TIFF; a vector file has no resolution
    dpi: int = 300
    font: str = "Arial"
    font_size: float = 8.0
    title: str = ""
    #: None: the content's own
    xlabel: str | None = None
    ylabel: str | None = None
    #: diameter in points
    point_size: float = 2.5
    opacity: float = 0.85
    #: open (left and bottom), box, none
    frame: str = "open"
    ticks: bool = True
    grid: bool = False
    #: right, bottom, inside, none
    legend: str = "right"
    legend_title: str | None = None
    legend_counts: bool = False
    #: category name → colour, where it differs from the default
    colors: dict[str, str] = field(default_factory=dict)
    #: category name → what to print instead
    labels: dict[str, str] = field(default_factory=dict)
    colormap: str = "viridis"
    #: the colours of the ramp when `colormap` is "custom", low to high
    colormap_colors: list[str] = field(default_factory=list)
    #: numbers in the heatmap's cells
    annotate: bool = True
    transparent: bool = False
    #: points as a bitmap inside a vector file: for maps of very many
    #: events, where a million circles make an SVG nobody can open
    rasterize: bool = False

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "FigureStyle":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in (d or {}).items() if k in known})

    def color(self, name: str, default: str) -> str:
        return self.colors.get(name) or default

    def label(self, name: str) -> str:
        return self.labels.get(name) or (name if name != "" else "—")


#: what a preset carries: the look, not this figure's own words
PRESET_FIELDS = [f.name for f in fields(FigureStyle)
                 if f.name not in ("title", "xlabel", "ylabel", "legend_title")]


# --- drawing -----------------------------------------------------------------------

def _figure(style: FigureStyle):
    from matplotlib.figure import Figure
    fig = Figure(figsize=(style.width_mm / 25.4, style.height_mm / 25.4),
                 layout="constrained")
    if style.transparent:
        fig.patch.set_alpha(0.0)
    return fig


def _rc(style: FigureStyle) -> dict:
    return {
        "font.family": [style.font, "DejaVu Sans"], "font.size": style.font_size,
        "axes.titlesize": style.font_size + 1, "axes.labelsize": style.font_size,
        "xtick.labelsize": style.font_size - 1, "ytick.labelsize": style.font_size - 1,
        "legend.fontsize": style.font_size - 1, "legend.title_fontsize": style.font_size - 1,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        # text stays text, and keeps its font
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.default": "regular", "axes.unicode_minus": False,
    }


def _frame(ax, style: FigureStyle) -> None:
    if style.transparent:
        ax.patch.set_alpha(0.0)
    if style.frame == "open":
        ax.spines[["top", "right"]].set_visible(False)
    elif style.frame == "none":
        ax.spines[:].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
    if not style.ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    if style.grid and style.frame != "none":
        ax.grid(True, linewidth=0.3, color="#d0d0d0")
        ax.set_axisbelow(True)


def _legend(fig, axes, handles, style: FigureStyle, title: str) -> None:
    if style.legend == "none" or not handles:
        return
    title = style.legend_title if style.legend_title is not None else title
    kw = dict(title=title or None, frameon=False, handletextpad=0.4,
              borderaxespad=0.2, labelspacing=0.35)
    if style.legend == "inside":
        axes[0].legend(handles=handles, loc="best", **kw)
    elif style.legend == "bottom":
        fig.legend(handles=handles, loc="outside lower center",
                   ncols=min(len(handles), 5), columnspacing=1.0, **kw)
    else:
        fig.legend(handles=handles, loc="outside right upper", **kw)


_SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def _tick_text(label: str) -> str:
    """`10³` as mathtext: few fonts carry every superscript digit, and a
    missing one prints as an empty box."""
    if label.startswith("10") and len(label) > 2 and label[2:].translate(_SUPERSCRIPT) != label[2:]:
        return "$10^{%s}$" % label[2:].translate(_SUPERSCRIPT)
    return label


def _dot(color: str, label: str):
    from matplotlib.lines import Line2D
    return Line2D([], [], marker="o", linestyle="", markersize=4.5,
                  markerfacecolor=color, markeredgewidth=0, label=label)


def _draw_map(fig, c: MapContent, style: FigureStyle) -> None:
    from matplotlib.colors import Normalize, to_rgba

    n = c.x.size
    order = np.random.default_rng(0).permutation(n)  # as on the screen
    groups = ([(None, np.ones(n, dtype=bool))] if c.split is None else
              [(name, c.split.codes == i) for i, name in enumerate(c.split.names)])
    # up to four panels read best side by side; more than that, a grid
    cols = c.columns or (len(groups) if len(groups) <= 4
                         else int(np.ceil(np.sqrt(len(groups)))))
    cols = max(1, min(cols, len(groups)))
    rows = -(-len(groups) // cols)
    grid = fig.subplots(rows, cols, sharex=True, sharey=True, squeeze=False)
    axes = [grid[k // cols][k % cols] for k in range(len(groups))]
    for extra in grid.ravel()[len(groups):]:
        extra.set_visible(False)

    mappable = None
    if isinstance(c.color, Categories):
        table = np.array([to_rgba(style.color(nm, col), style.opacity)
                          for nm, col in zip(c.color.names, c.color.colors)])
        rgba = table[c.color.codes]
    elif isinstance(c.color, Values):
        lo, hi = c.color.limits()
        norm = Normalize(lo, hi, clip=True)
        cmap = colormap(style.colormap, style.colormap_colors)
        rgba = cmap(norm(np.nan_to_num(c.color.values, nan=lo)))
        rgba[:, 3] = style.opacity
        from matplotlib.cm import ScalarMappable
        mappable = ScalarMappable(norm, cmap)
    else:
        rgba = np.tile(to_rgba(style.color("", "#234b8c"), style.opacity), (n, 1))

    area = style.point_size ** 2
    for ax, (name, member) in zip(axes, groups):
        if name is not None:
            rest = order[~member[order]]
            ax.scatter(c.x[rest], c.y[rest], s=area, c="#d9dbe0", linewidths=0,
                       rasterized=style.rasterize)
            count = f" ({int(member.sum()):,d})" if style.legend_counts else ""
            ax.set_title(style.label(name) + count)
        sel = order[member[order]]
        ax.scatter(c.x[sel], c.y[sel], s=area, c=rgba[sel], linewidths=0,
                   rasterized=style.rasterize)
        if c.xticks:
            ax.set_xticks([p for p, _ in c.xticks], [_tick_text(t) for _, t in c.xticks])
        if c.yticks:
            ax.set_yticks([p for p, _ in c.yticks], [_tick_text(t) for _, t in c.yticks])
        _frame(ax, style)
    xlabel = c.xlabel if style.xlabel is None else style.xlabel
    ylabel = c.ylabel if style.ylabel is None else style.ylabel
    if style.frame != "none":
        for k, ax in enumerate(axes):
            if k // cols == rows - 1 or k + cols >= len(groups):
                ax.set_xlabel(xlabel)
            if k % cols == 0:
                ax.set_ylabel(ylabel)

    if isinstance(c.color, Categories):
        counts = c.color.counts()
        handles = [_dot(style.color(nm, col),
                        style.label(nm) + (f" ({int(k):,d})" if style.legend_counts else ""))
                   for nm, col, k in zip(c.color.names, c.color.colors, counts)]
        _legend(fig, axes, handles, style, c.color.title)
    elif mappable is not None and style.legend != "none":
        bar = fig.colorbar(mappable, ax=axes, shrink=0.7, aspect=22, pad=0.02,
                           location="bottom" if style.legend == "bottom" else "right")
        bar.set_label(style.legend_title if style.legend_title is not None
                      else c.color.title)
        if c.color.ticks:
            inside = [(p, _tick_text(t)) for p, t in c.color.ticks
                      if mappable.norm.vmin <= p <= mappable.norm.vmax]
            if len(inside) >= 2:
                bar.set_ticks([p for p, _ in inside], labels=[t for _, t in inside])
        bar.outline.set_linewidth(0.5)


def _draw_heatmap(fig, c: HeatmapContent, style: FigureStyle) -> None:
    t = c.table
    shade = cluster_stats.scaled(t, c.per_column)
    ax = fig.subplots()
    cmap = colormap(style.colormap, style.colormap_colors)
    im = ax.imshow(shade, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(t.columns)),
                  [style.label(n.split()[0] if n else n) for n in t.columns],
                  rotation=45 if len(t.columns) > 6 else 0,
                  ha="right" if len(t.columns) > 6 else "center")
    ax.set_yticks(range(len(t.rows)), [style.label(r) for r in t.rows])
    ax.set_ylabel("Cluster" if style.ylabel is None else style.ylabel)
    ax.set_xlabel("" if style.xlabel is None else style.xlabel)
    ax.tick_params(length=0)
    ax.spines[:].set_visible(style.frame == "box")
    if style.annotate:
        for r in range(shade.shape[0]):
            for k in range(shade.shape[1]):
                red, green, blue, _ = cmap(float(np.clip(shade[r, k], 0, 1)))
                dark = 0.299 * red + 0.587 * green + 0.114 * blue < 0.55
                ax.text(k, r, f"{t.values[r, k]:.0f}", ha="center", va="center",
                        fontsize=style.font_size - 1.5,
                        color="white" if dark else "black")
    if style.legend != "none":
        bar = fig.colorbar(im, ax=ax, shrink=0.7, aspect=22, pad=0.02, ticks=[0, 1],
                           location="bottom" if style.legend == "bottom" else "right")
        bar.ax.set_yticklabels(["min", "max"]) if style.legend != "bottom" else \
            bar.ax.set_xticklabels(["min", "max"])
        default = "Per marker" if c.per_column else "Whole table"
        bar.set_label(style.legend_title if style.legend_title is not None else default)
        bar.outline.set_linewidth(0.5)


def _draw_frequencies(fig, c: FrequencyContent, style: FigureStyle) -> None:
    panels = cluster_stats.ordered(c.facets) if c.facets else [None]
    axes = fig.subplots(1, len(panels), sharey=True, squeeze=False)[0]
    handles = []
    for k, (ax, panel) in enumerate(zip(axes, panels)):
        keep = [i for i in range(len(c.groups)) if panel is None or c.facets[i] == panel]
        handles = _frequency_panel(ax, c, style, keep, first=k == 0) or handles
        if panel is not None:
            ax.set_title(style.label(panel))
    _legend(fig, list(axes), handles, style, "")


def _frequency_panel(ax, c: FrequencyContent, style: FigureStyle, keep: list[int],
                     first: bool) -> list:
    """One panel of bars: the samples in `keep`, by group."""
    from matplotlib.patches import Patch
    names = cluster_stats.ordered(c.groups)      # every group in every panel:
    colors = categorical(len(names))             # a colour means one thing
    n_clusters = len(c.freq.rows)
    width = 0.8 / max(len(names), 1)
    rng = np.random.default_rng(0)
    handles = []
    for j, (g, default) in enumerate(zip(names, colors)):
        color = style.color(g, default)
        cols = [k for k in keep if c.groups[k] == g]
        vals = c.freq.values[:, cols]
        vals = vals[:, ~np.isnan(vals).all(axis=0)]
        if vals.shape[1] == 0:
            handles.append(Patch(facecolor=color, alpha=0.6, edgecolor=color,
                                 label=style.label(g)))
            continue
        x = np.arange(n_clusters) - 0.4 + width * (j + 0.5)
        mean = vals.mean(axis=1)
        ax.bar(x, mean, width * 0.9, color=color, alpha=0.6, edgecolor=color, linewidth=0.6)
        if vals.shape[1] > 1:
            ax.errorbar(x, mean, yerr=vals.std(axis=1, ddof=1), fmt="none",
                        ecolor="black", elinewidth=0.6, capsize=1.5, capthick=0.6)
        jitter = (rng.random(vals.shape) - 0.5) * width * 0.5
        ax.scatter((x[:, None] + jitter).ravel(), vals.ravel(),
                   s=max(style.point_size, 1.5) ** 2, facecolors="white",
                   edgecolors="black", linewidths=0.5, zorder=3)
        suffix = f" (n={vals.shape[1]})" if style.legend_counts else ""
        handles.append(Patch(facecolor=color, alpha=0.6, edgecolor=color,
                             label=style.label(g) + suffix))
    ax.set_xticks(range(n_clusters), [style.label(r) for r in c.freq.rows])
    ax.set_xlabel("Cluster" if style.xlabel is None else style.xlabel)
    if first:
        ax.set_ylabel("% of the sample's events" if style.ylabel is None else style.ylabel)
    ax.set_xlim(-0.6, n_clusters - 0.4)
    ax.set_ylim(bottom=0)
    _frame(ax, style)
    return handles


def _draw_tree(fig, c: TreeContent, style: FigureStyle) -> None:
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    t = c.tree
    ax = fig.subplots()
    pos, sizes, cluster = t["pos"], t["sizes"], t["cluster"]
    used = np.flatnonzero(sizes > 0)
    if len(t["edges"]):
        ax.add_collection(LineCollection(pos[t["edges"]], colors="#b8bcc4",
                                         linewidths=0.7, zorder=1))
    area = (style.point_size * (1.2 + 3.5 * np.sqrt(sizes[used] / sizes[used].max()))) ** 2
    if c.color_by in t["markers"]:
        m = t["medians"][used, t["markers"].index(c.color_by)]
        norm = Normalize(np.nanmin(m), np.nanmax(m))
        cmap = colormap(style.colormap, style.colormap_colors)
        ax.scatter(pos[used, 0], pos[used, 1], s=area, c=m, cmap=cmap,
                   norm=norm, edgecolors="white", linewidths=0.4, zorder=2)
        if style.legend != "none":
            bar = fig.colorbar(ScalarMappable(norm, cmap), ax=ax,
                               shrink=0.6, aspect=22, pad=0.02,
                               location="bottom" if style.legend == "bottom" else "right")
            bar.set_label(style.legend_title if style.legend_title is not None
                          else c.color_by.split()[0])
            bar.outline.set_linewidth(0.5)
    else:
        named = dict(category_names(c))
        colors = [style.color(str(int(k)), named[str(int(k))]) for k in cluster[used]]
        ax.scatter(pos[used, 0], pos[used, 1], s=area, c=colors,
                   edgecolors="white", linewidths=0.4, zorder=2)
        handles = [_dot(style.color(nm, col), style.label(nm)) for nm, col in named.items()]
        _legend(fig, [ax], handles, style, "Cluster")
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.08)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(style.frame == "box")
    if style.transparent:
        ax.patch.set_alpha(0.0)


def _draw_overlay(fig, c: OverlayContent, style: FigureStyle) -> None:
    from matplotlib.lines import Line2D
    ax = fig.subplots()
    colors = dict(category_names(c))
    colour = {g: style.color(g, default) for g, default in colors.items()}
    if c.is_histogram:
        top = max((float(s.y.max()) for s in c.series if s.y.size), default=1.0) or 1.0
        step = STAGGER * top if c.offset else 0.0
        for k, s in enumerate(c.series):
            base = k * step
            ax.fill_between(s.x, base, s.y + base, color=colour[s.group],
                            alpha=0.18 if not c.offset else 0.35, linewidth=0,
                            zorder=2 * (len(c.series) - k))
            ax.plot(s.x, s.y + base, color=colour[s.group], linewidth=0.9,
                    zorder=2 * (len(c.series) - k) + 1)
        if c.offset:
            ax.set_yticks([k * step for k in range(len(c.series))],
                          [style.label(s.label) for s in c.series])
            ax.tick_params(axis="y", length=0)
            ax.set_ylim(0, (len(c.series) - 1) * step + top * 1.05)
        else:
            ax.set_ylim(0, top * 1.05)
            ax.set_ylabel(c.ylabel if style.ylabel is None else style.ylabel)
    else:
        widths = {0.90: 0.5, 0.75: 0.7, 0.50: 0.9, 0.25: 1.1, 0.10: 1.3}
        for s in c.series:
            for share, lines in s.levels or []:
                for line in lines:
                    ax.plot(line[:, 0], line[:, 1], color=colour[s.group],
                            linewidth=widths.get(round(share, 2), 0.8))
            if s.outliers is not None and len(s.outliers):
                ax.scatter(s.outliers[:, 0], s.outliers[:, 1], s=(style.point_size * 0.5) ** 2,
                           color=colour[s.group], linewidths=0, alpha=0.6,
                           rasterized=style.rasterize)
        if c.yticks:
            ax.set_yticks([p for p, _ in c.yticks], [_tick_text(t) for _, t in c.yticks])
        if c.yrange:
            ax.set_ylim(*c.yrange)
        ax.set_ylabel(c.ylabel if style.ylabel is None else style.ylabel)
    if c.xticks:
        ax.set_xticks([p for p, _ in c.xticks], [_tick_text(t) for _, t in c.xticks])
    if c.xrange:
        ax.set_xlim(*c.xrange)
    ax.set_xlabel(c.xlabel if style.xlabel is None else style.xlabel)
    _frame(ax, style)
    handles = [Line2D([], [], color=colour[g], linewidth=1.6, label=style.label(g))
               for g in colors]
    # staggered and one colour per sample: the axis already names every curve
    if not (c.is_histogram and c.offset and len(colors) == len(c.series)):
        _legend(fig, [ax], handles, style, "")


_DRAW = {OverlayContent: _draw_overlay, MapContent: _draw_map, HeatmapContent: _draw_heatmap,
         FrequencyContent: _draw_frequencies, TreeContent: _draw_tree}


class _drawn:
    """`with _drawn(content, style) as fig:` — the figure, built and still
    under the style's settings, which matplotlib reads again when it
    writes (the SVG font handling among them)."""

    def __init__(self, content, style: FigureStyle):
        self.content, self.style = content, style

    def __enter__(self):
        import matplotlib
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        self._rc = matplotlib.rc_context(_rc(self.style))
        self._rc.__enter__()
        fig = _figure(self.style)
        _DRAW[type(self.content)](fig, self.content, self.style)
        if self.style.title:
            fig.suptitle(self.style.title, fontsize=self.style.font_size + 1.5)
        FigureCanvasAgg(fig)
        return fig

    def __exit__(self, *exc):
        return self._rc.__exit__(*exc)


def save(content, style: FigureStyle, path: str) -> None:
    """The format is the file's extension: .svg, .pdf, .png, .tif(f)."""
    ext = path.rsplit(".", 1)[-1].lower()
    fmt = {"tif": "tiff"}.get(ext, ext)
    if fmt not in FORMATS:
        raise ValueError(f"cannot write .{ext}: use " + ", ".join("." + f for f in FORMATS))
    kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
    with _drawn(content, style) as fig:
        fig.savefig(path, format=fmt, dpi=style.dpi, transparent=style.transparent, **kw)


def preview_png(content, style: FigureStyle, dpi: int = 110) -> bytes:
    """The figure as it will be exported, small: drawn by the code that
    writes the file, so the preview cannot disagree with it."""
    import io
    buf = io.BytesIO()
    with _drawn(content, style) as fig:
        fig.savefig(buf, format="png", dpi=dpi, transparent=False)
    return buf.getvalue()


def fonts() -> list[str]:
    """The font families matplotlib can actually draw with here."""
    from matplotlib import font_manager
    return sorted(set(font_manager.get_font_names()), key=str.lower)
