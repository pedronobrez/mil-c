"""
Several panels on one page, as a paper figure is, and a report: the pages
one after the other in one PDF, behind a first page that says what the
analysis was.

A panel is a *content* of `figures.py` and is drawn by the function that
draws it alone, into its own cell of the page (a matplotlib sub-figure),
so a panel cannot look different from the single figure it came from.
The page has one `FigureStyle`: its type is every panel's, which is what
makes six plots one figure. What a panel says of itself — its labels, its
legend, its colours — it keeps.

No Qt here, and matplotlib is imported only when a page is drawn.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass, field, replace

from . import figures
from .figures import FigureStyle

#: what belongs to the page and is never a panel's own: one type and one
#: resolution across the panels (and the size, which is the page's)
PAGE_FIELDS = ("width_mm", "height_mm", "dpi", "font", "font_size", "transparent", "title")

#: how much taller than the type a panel's letter is, in points
LETTER_EXTRA = 4.0

#: the page of a report that is text: A4, upright
REPORT_PAGE_MM = (210.0, 297.0)
REPORT_MARGIN_MM = 20.0
_HEADINGS = ("Sampling", "Batch normalisation", "Runs")


@dataclass
class Panel:
    content: object
    title: str = ""
    #: `FigureStyle` fields where this panel differs from the page
    style: dict = field(default_factory=dict)
    #: how many columns and rows of the grid it takes
    colspan: int = 1
    rowspan: int = 1
    #: where it came from, to ask again when the analysis has changed:
    #: anything with `content()`
    source: object = None
    #: what a list of panels calls it: `map`, `heatmap`…
    name: str = ""
    #: the source had nothing to show at the last refresh: this is the
    #: content of before
    stale: bool = False


def own_style(style: FigureStyle, page: FigureStyle) -> dict:
    """What a panel drawn with `style` must keep of it to look the same on
    a page whose style is `page`."""
    theirs, ours = style.to_dict(), page.to_dict()
    return {k: v for k, v in theirs.items() if k not in PAGE_FIELDS and v != ours[k]}


@dataclass
class Plate:
    panels: list[Panel] = field(default_factory=list)
    #: cells per row of the grid
    columns: int = 2
    width_mm: float = 180.0
    #: None: three quarters of a cell's width for each row
    height_mm: float | None = None
    #: the type, the resolution and the default look of every panel; its
    #: own size and title are not used
    style: FigureStyle = field(default_factory=FigureStyle)
    #: A, B, C… at the top left of each panel, bold, `LETTER_EXTRA` points
    #: above the style's font size
    letters: bool = True

    def add(self, content, style: FigureStyle | None = None, title: str | None = None,
            source=None, name: str = "") -> Panel:
        """A panel as `style` draws it alone; the style's title becomes the
        panel's unless another is given."""
        style = style or self.style
        panel = Panel(content, style.title if title is None else title,
                      own_style(style, self.style), source=source, name=name)
        self.panels.append(panel)
        return panel

    def move(self, index: int, by: int) -> int:
        """Where the panel is afterwards: the ends are walls."""
        to = min(max(index + by, 0), len(self.panels) - 1)
        self.panels.insert(to, self.panels.pop(index))
        return to

    def style_of(self, panel: Panel) -> FigureStyle:
        own = {k: v for k, v in panel.style.items() if k not in PAGE_FIELDS}
        for key in ("colors", "labels"):   # a name's colour on the page, unless the panel says
            own[key] = {**getattr(self.style, key), **own.get(key, {})}
        return replace(self.style, **own, title="")

    def cells(self) -> list[tuple[int, int, int, int]]:
        """(row, column, rows, columns) of each panel: each in the first
        free place after the one before that holds it, so the letters
        read across and down whatever the panels span."""
        taken: set[tuple[int, int]] = set()
        out = []
        k = 0
        for p in self.panels:
            cs = min(max(p.colspan, 1), self.columns)
            rs = max(p.rowspan, 1)
            while True:
                r, c = divmod(k, self.columns)
                block = {(r + i, c + j) for i in range(rs) for j in range(cs)}
                if c + cs <= self.columns and not block & taken:
                    break
                k += 1
            taken |= block
            out.append((r, c, rs, cs))
        return out

    def rows(self) -> int:
        return max((r + rs for r, _, rs, _ in self.cells()), default=1)

    def size_mm(self) -> tuple[float, float]:
        if self.height_mm is not None:
            return self.width_mm, self.height_mm
        return self.width_mm, round(self.rows() * self.width_mm / self.columns * 0.75, 1)

    def refresh(self) -> list[Panel]:
        """Ask every panel's source what it shows now. A source with
        nothing to show leaves the panel as it was, marked stale; those
        panels are returned."""
        for p in self.panels:
            if p.source is None:
                continue
            content = p.source.content()
            p.stale = content is None
            if content is not None:
                p.content = content
        return [p for p in self.panels if p.stale]


def letter(index: int) -> str:
    """A to Z, then AA, AB…"""
    text = ""
    index += 1
    while index:
        index, rest = divmod(index - 1, 26)
        text = chr(65 + rest) + text
    return text


# --- drawing -----------------------------------------------------------------------

def _draw(fig, plate: Plate) -> None:
    from matplotlib.transforms import offset_copy
    grid = fig.add_gridspec(plate.rows(), plate.columns)
    pad = fig.get_layout_engine().get()["h_pad"]
    size = plate.style.font_size
    for k, (panel, (r, c, rs, cs)) in enumerate(zip(plate.panels, plate.cells())):
        sub = fig.add_subfigure(grid[r:r + rs, c:c + cs])
        figures._DRAW[type(panel.content)](sub, panel.content, plate.style_of(panel))
        if not plate.letters:
            if panel.title:
                sub.suptitle(panel.title, fontsize=size + 1.5)
            continue
        # the letter is the sub-figure's title, set to the left: the layout
        # then keeps a band for it above the plot, and the panel's own
        # title is written in that band
        inches = sub.bbox.width / fig.dpi
        sub.suptitle(letter(k), x=pad / inches, ha="left", fontweight="bold",
                     fontsize=size + LETTER_EXTRA)
        if panel.title:
            top = offset_copy(sub.transSubfigure, fig, y=-pad, units="inches")
            sub.text(0.5, 1.0, panel.title, ha="center", va="top", fontsize=size + 1.5,
                     transform=top)


class _drawn:
    """`with _drawn(page) as fig:` for a plate, or for a single figure
    given as (content, style): built, and still under the style's
    settings, which matplotlib reads again when it writes."""

    def __init__(self, page):
        self.page = page

    def __enter__(self):
        if not isinstance(self.page, Plate):
            self._inner = figures._drawn(*self.page)
            return self._inner.__enter__()
        import matplotlib
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        plate = self.page
        if not plate.panels:
            raise ValueError("the plate has no panels")
        self._inner = matplotlib.rc_context(figures._rc(plate.style))
        self._inner.__enter__()
        try:
            fig = figures._figure(replace(plate.style, width_mm=plate.size_mm()[0],
                                          height_mm=plate.size_mm()[1]))
            _draw(fig, plate)
            FigureCanvasAgg(fig)
        except BaseException:
            self._inner.__exit__(None, None, None)
            raise
        return fig

    def __exit__(self, *exc):
        return self._inner.__exit__(*exc)


def save(plate: Plate, path: str) -> None:
    """The format is the file's extension, as for a single figure."""
    ext = path.rsplit(".", 1)[-1].lower()
    fmt = {"tif": "tiff"}.get(ext, ext)
    if fmt not in figures.FORMATS:
        raise ValueError(f"cannot write .{ext}: use "
                         + ", ".join("." + f for f in figures.FORMATS))
    kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
    with _drawn(plate) as fig:
        fig.savefig(path, format=fmt, dpi=plate.style.dpi,
                    transparent=plate.style.transparent, **kw)


def preview_png(plate: Plate, dpi: int = 110) -> bytes:
    """The page as it will be exported, small, by the code that writes it."""
    import io
    buf = io.BytesIO()
    with _drawn(plate) as fig:
        fig.savefig(buf, format="png", dpi=dpi, transparent=False)
    return buf.getvalue()


# --- the report --------------------------------------------------------------------

_MODES = {"equal_per_sample": "the same number of events from every sample",
          "equal_per_group": "the same total from every group",
          "ceiling": "at most a fixed number from each sample"}


def _run_line(name: str, run: dict, dataset) -> str:
    known = set(dataset.channels)
    on = ", ".join(dataset.label(c) if c in known else c for c in run.get("channels", []))
    scale = run.get("transform", "")
    if "cofactor" in run:
        scale += f" (cofactor {run['cofactor']:g})"
    parts = [f"on {on}" if on else "", f"scale: {scale}" if scale else ""]
    if name == "FlowSOM":
        parts.append(f"{run.get('n_clusters')} clusters on a "
                     f"{run.get('grid')} × {run.get('grid')} grid")
    else:
        parts += [f"{k} {run[k]:g}" for k in ("n_neighbors", "min_dist", "perplexity")
                  if k in run]
    if "seed" in run:
        parts.append(f"seed {run['seed']}")
    line = f"{name}: " + "; ".join(p for p in parts if p) + "."
    if "neighbours_kept" in run:
        line += (f" Keeps {run['neighbours_kept']:.1%} of each event's 15 nearest "
                 "neighbours.")
    return line


def _normalisation_lines(done: dict) -> list[str]:
    column, value = done.get("reference", ["", ""])
    lines = [f"Batches by the column {done.get('batch')!r}; reference samples: "
             f"{column} = {value}.",
             "Markers: " + ", ".join(done.get("markers", [])) + ".",
             f"FlowSOM clusters within which it was done: {done.get('clusters_used', 1)}."]
    for batch in done.get("batches", []):
        names = ", ".join(done.get("reference_samples", {}).get(batch, []))
        events = done.get("reference_events", {}).get(batch, 0)
        lines.append(f"  Batch {batch or '—'}: {events:,d} reference events ({names})")
    before, after = done.get("distance_before"), done.get("distance_after")
    if before and after:
        lines.append("Distance between the batches' references, mean over markers, in "
                     f"scale units, on events held out: {sum(before) / len(before):.4f} "
                     f"before, {sum(after) / len(after):.4f} after.")
    for batch, cluster, events in done.get("fallbacks", []):
        lines.append(f"  Batch {batch}, cluster {cluster}: {events} reference events; "
                     "the whole batch's curve was used.")
    if done.get("not_normalised"):
        lines.append("Not normalised (no reference in their batch): "
                     + ", ".join(done["not_normalised"]) + ".")
    return lines


def summary(dataset) -> list[str]:
    """What the analysis was, as the dataset records it, one line each;
    an empty string is a blank line. Only what is written there: nothing
    is recomputed and nothing is interpreted."""
    p = dataset.provenance
    lines = [f"{len(dataset):,d} events from {len(dataset.samples)} sample(s)."]
    if p.get("source"):
        lines.append(f"Read from an already concatenated file: {p['source']}.")
    if "population" in p:
        lines.append(f"Population: {p['population'] or 'all events'}.")
    if "mode" in p:
        lines += ["", _HEADINGS[0],
                  f"{_MODES.get(p['mode'], p['mode']).capitalize()}; seed {p.get('seed')}."]
        if p.get("limited_by"):
            lines.append(f"The number was set by {p['limited_by']}.")
        available = p.get("available", {})
        for name, taken in p.get("take", {}).items():
            lines.append(f"  {name}: {taken:,d} of {available.get(name, taken):,d}")
        for name, why in p.get("excluded", {}).items():
            lines.append(f"  {name}: left out ({why})")
    if p.get("normalization"):
        lines += ["", _HEADINGS[1]] + _normalisation_lines(p["normalization"])
    runs = [name for name in ("UMAP", "tSNE", "FlowSOM") if name in dataset.runs]
    runs += [name for name in dataset.runs if name not in runs]
    if runs:
        lines += ["", _HEADINGS[2]] + [_run_line(name, dataset.runs[name], dataset)
                                 for name in runs]
    return lines


def _text_pages(lines: list[str], title: str, style: FigureStyle):
    """The lines as A4 pages of plain text, as many as they need."""
    from matplotlib.figure import Figure
    width, height = REPORT_PAGE_MM
    size = style.font_size + 1
    step = size * 1.5 / 72 * 25.4 / height            # one line, in page heights
    # a proportional font: about half its size per character, on average
    chars = int((width - 2 * REPORT_MARGIN_MM) / 25.4 * 72 / (size * 0.5))
    wrapped = []
    for line in lines:
        indent = " " * (len(line) - len(line.lstrip()) + 2)
        wrapped += textwrap.wrap(line, chars, subsequent_indent=indent) or [""]
    x, top = REPORT_MARGIN_MM / width, 1 - REPORT_MARGIN_MM / height
    room = int((top - REPORT_MARGIN_MM / height) / step) - 3
    for start in range(0, max(len(wrapped), 1), room):
        fig = Figure(figsize=(width / 25.4, height / 25.4))
        fig.text(x, top, title, fontsize=size + 5, fontweight="bold", va="top")
        for k, line in enumerate(wrapped[start:start + room]):
            indent = len(line) - len(line.lstrip())
            fig.text(x + 0.01 * indent, top - (k + 3) * step, line.lstrip(),
                     fontsize=size, va="top",
                     fontweight="bold" if line in _HEADINGS else "normal")
        yield fig


def report(path: str, pages: list, dataset=None, title: str = "MIL-C analysis") -> int:
    """A PDF of several pages: each of `pages` — a `Plate`, or one figure
    as (content, style) — on a page of its own size, behind what
    `summary` says of `dataset` when one is given. The number of pages
    written."""
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages
    first = pages[0] if pages else None
    style = (first.style if isinstance(first, Plate) else first[1]) if first else FigureStyle()
    count = 0
    # the fonts are embedded when the file is closed: the settings that keep
    # text as text must still hold then
    with matplotlib.rc_context(figures._rc(style)), PdfPages(path) as pdf:
        if dataset is not None:
            for fig in _text_pages(summary(dataset), title, style):
                pdf.savefig(fig)
                count += 1
        for page in pages:
            with _drawn(page) as fig:
                pdf.savefig(fig)
                count += 1
    return count
