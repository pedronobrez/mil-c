"""
The gates, made ready to be written for another program: what a FlowJo
workspace (`wsp_write.py`) and a Gating-ML file (`gatingml.py`) have in
common.

Both formats hold rectangles, polygons, ellipses, quadrants and
combinations, and nothing else. So a tree is first rewritten with those
alone (`prepare`): a time gate becomes a range on Time — several ranges
and their union where it keeps several stretches —, a population that is
none of several becomes the negation of their union. What cannot be said
that way (a cluster of a map, a kind of gate not known here) is left out
with everything under it, and listed in the `Report`. Nothing here raises
for one gate: an export loses that gate, not the others.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .compensation import Compensation, strip_prefix
from .gating import (AND, NOT, OR, BooleanGate, ClusterGate, EllipseGate, GateNode,
                     PolygonGate, QuadrantGate, RectangleGate, TimeGate)
from .transforms import Biex, Linear, Logicle, same, to_spec

WSP, GATINGML = "wsp", "gatingml"
FORMATS = {WSP: "FlowJo workspace", GATINGML: "Gating-ML 2.0"}

#: how far a side redrawn on another scale may be from the one drawn, as a
#: fraction of the polygon's own size
TOLERANCE = 1e-5
#: the sides of the polygon an ellipse becomes where it cannot stay one
ELLIPSE_SIDES = 512


@dataclass
class Skipped:
    """A population that was not written."""
    sample: str
    path: str
    #: the gate's class here: `ClusterGate`, `TimeGate`…
    kind: str
    reason: str
    #: how many populations under it went with it
    below: int = 0


@dataclass
class Note:
    """A population written in another form than it has here."""
    sample: str
    path: str
    text: str


@dataclass
class Report:
    format: str = ""
    files: list[str] = field(default_factory=list)
    samples: int = 0
    #: populations of the trees that were written, helpers not counted
    written: int = 0
    skipped: list[Skipped] = field(default_factory=list)
    notes: list[Note] = field(default_factory=list)

    def skip(self, sample: str, node: GateNode, reason: str) -> None:
        self.skipped.append(Skipped(sample, source_of(node), type(node.gate).__name__, reason,
                                    sum(1 for _ in node.walk()) - 1))

    def note(self, sample: str, path: str, text: str) -> None:
        if not any((n.sample, n.path, n.text) == (sample, path, text) for n in self.notes):
            self.notes.append(Note(sample, path, text))

    def summary(self) -> str:
        """One line: what was written, and how much was not."""
        what = (os.path.basename(self.files[0]) if len(self.files) == 1
                else f"{len(self.files)} files")
        text = (f"Wrote {what}: {self.written} population(s) of {self.samples} sample(s)")
        if self.skipped:
            lost = sum(1 + s.below for s in self.skipped)
            text += f"; {lost} not written"
        return text

    def text(self) -> str:
        """Everything a reader should know, in lines: the files, what was
        left out and why, what was written in another form."""
        lines = [self.summary() + "."]
        if len(self.files) > 1:
            lines += ["  " + f for f in self.files]
        for title, rows in (("Not written", [(s.sample, s.path, s.reason
                                              + (f" (and {s.below} under it)" if s.below else ""))
                                             for s in self.skipped]),
                            ("Written in another form", [(n.sample, n.path, n.text)
                                                         for n in self.notes])):
            if not rows:
                continue
            lines += ["", title + ":"]
            # the same thing said of every sample is said once
            grouped: dict[tuple, list[str]] = {}
            for sample, path, why in rows:
                grouped.setdefault((path, why), []).append(sample)
            for (path, why), samples in grouped.items():
                who = (samples[0] if len(samples) == 1 else f"{len(samples)} samples")
                lines.append(f"  {path or 'every gate'} ({who}): {why}" if who
                             else f"  {path or 'every gate'}: {why}")
        return "\n".join(lines)


@dataclass
class ExportSample:
    """One file and what its gates are applied to."""
    #: the file's name, with its extension: what FlowJo calls the sample
    name: str
    path: str
    events: int
    channels: list[str]
    keywords: dict[str, str] = field(default_factory=dict)
    compensation: Compensation | None = None
    #: `file`, `workspace`, `edited`, `shared` or ""
    matrix_source: str = ""
    #: channel → the scale of its axis; a channel not here, or None, is
    #: drawn as stored
    scales: dict[str, object] = field(default_factory=dict)
    gates: GateNode = field(default_factory=lambda: GateNode(""))
    #: gives the opened `Sample`, for the counts; None writes none
    load: Callable | None = None

    def dimension(self, name: str) -> tuple[str, bool]:
        """(the file's channel, is it compensated) for a gate's dimension,
        as `Sample.column` reads it."""
        if name in self.channels:
            return name, False
        base, wants = strip_prefix(name)
        if base not in self.channels:
            raise Unwritable(f"the file has no channel {name}")
        return base, bool(wants and self.compensation is not None
                          and base in self.compensation.channels)

    def range_of(self, channel: str) -> float:
        """`$PnR`: the top of a channel as the file declares it."""
        try:
            return float(self.keywords[f"$P{self.channels.index(channel) + 1}R"])
        except (KeyError, ValueError):
            return 262144.0

    @property
    def time(self) -> str | None:
        return next((c for c in self.channels if c.lower() == "time"), None)

    @property
    def timestep(self) -> float:
        """Seconds in one unit of the file's Time, 1 where it does not say."""
        try:
            step = float(self.keywords.get("$TIMESTEP", "") or 1.0)
        except ValueError:
            step = 1.0
        return step if step > 0 and np.isfinite(step) else 1.0


@dataclass
class Prepared:
    sample: ExportSample
    #: the tree with only gates both formats hold; see `prepare`
    tree: GateNode
    #: path in `tree` → events, where the sample could be opened
    counts: dict[str, int] | None = None


class Unwritable(Exception):
    """A gate that is left out; the message is the reason."""


def source_of(node: GateNode) -> str:
    """The path the population has in MIL-C; for one added by `prepare`,
    its own."""
    return getattr(node, "source", None) or node.path


def is_helper(node: GateNode) -> bool:
    return getattr(node, "helper", False)


def _is(gate, kind) -> bool:
    """`gate` is a `kind`, or a kind of it that tests events the same way
    — not one that only shares its fields."""
    return isinstance(gate, kind) and getattr(type(gate), "contains", None) is getattr(
        kind, "contains", None) and getattr(type(gate), "combine", None) is getattr(
        kind, "combine", None)


# --- a tree both formats can hold --------------------------------------------------

def time_ranges(gate: TimeGate, double: bool = False) -> list[tuple[float | None, float | None]]:
    """A time gate as ranges `lo <= t < hi`, None for an open side: the
    stretches it keeps. Its intervals include both ends and a range does
    not include its top, so the top is the next value a file can hold.
    A file stores Time in single precision (`double` where it does not),
    and that is the precision events are compared in: the bounds are
    numbers of that precision, the same in whatever unit they are then
    written."""
    kind = np.float64 if double else np.float32
    spans = sorted((float(kind(lo)), float(kind(hi))) for lo, hi in gate.intervals if lo <= hi)
    merged: list[list[float]] = []
    for lo, hi in spans:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])

    def after(hi: float) -> float:
        return float(np.nextafter(kind(hi), kind(np.inf)))

    if gate.inside:
        return [(lo, after(hi)) for lo, hi in merged]
    edges = [None] + [v for lo, hi in merged for v in (lo, after(hi))] + [None]
    return list(zip(edges[::2], edges[1::2]))


def _unique(name: str, taken: set[str]) -> str:
    out, k = name, 2
    while out in taken:
        out, k = f"{name} ({k})", k + 1
    taken.add(out)
    return out


def _checked(gate):
    """The gate, if its numbers can be written at all."""
    def finite(values) -> bool:
        return bool(np.all(np.isfinite(np.asarray(values, dtype=np.float64))))

    if _is(gate, RectangleGate):
        if not 1 <= len(gate.dims) <= 2 or len(gate.mins) != len(gate.dims) \
                or len(gate.maxs) != len(gate.dims):
            raise Unwritable("a rectangle is on one or two parameters")
        if not finite([v for v in (*gate.mins, *gate.maxs) if v is not None]):
            raise Unwritable("one of its bounds is not a number")
    elif _is(gate, PolygonGate):
        v = np.asarray(gate.vertices, dtype=np.float64)
        if len(gate.dims) != 2 or v.ndim != 2 or v.shape[1] != 2 or len(v) < 3:
            raise Unwritable("a polygon has three vertices or more, on two parameters")
        if not finite(v):
            raise Unwritable("one of its vertices is not a number")
    elif _is(gate, EllipseGate):
        if len(gate.dims) != 2 or not finite([*gate.center, *gate.radii, gate.angle]) \
                or min(gate.radii) <= 0:
            raise Unwritable("an ellipse has a centre and two half-axes above zero")
    elif _is(gate, QuadrantGate):
        if len(gate.dims) != 2 or not finite([gate.x, gate.y]):
            raise Unwritable("a quadrant is cut at two numbers, on two parameters")
    return gate


def prepare(sample: ExportSample, report: Report) -> GateNode:
    """The sample's tree with only rectangles, polygons, ellipses,
    quadrants and combinations in it. A population keeps its name and its
    place; each carries `source`, its path in MIL-C. One added to say
    another (`helper`) has none of its own."""
    root = GateNode(sample.gates.name)

    def expand(node: GateNode, name: str, taken: set[str]) -> list[GateNode]:
        g = node.gate
        if isinstance(g, ClusterGate):
            raise Unwritable("a cluster is the events nearest some nodes of a map: "
                             "no other program has the map")
        if _is(g, TimeGate):
            if not g.dims:
                raise Unwritable("it names no Time parameter")
            ranges = time_ranges(g, sample.keywords.get("$DATATYPE", "").strip().upper() == "D")
            dim = g.dims[0]
            if not ranges:               # keeps nothing: an empty range says so
                ranges = [(0.0, 0.0)]
            if len(ranges) == 1:
                return [GateNode(name, RectangleGate([dim], [ranges[0][0]], [ranges[0][1]]))]
            parts = [GateNode(_unique(f"{name} [{k}]", taken),
                              RectangleGate([dim], [lo], [hi]))
                     for k, (lo, hi) in enumerate(ranges, 1)]
            for p in parts:
                p.helper = True
            report.note(sample.name, node.path,
                        f"the {len(ranges)} stretches of time it keeps are written as "
                        f"{len(ranges)} ranges beside it, and it as their union")
            return parts + [GateNode(name, BooleanGate(OR, [p.name for p in parts]))]
        if _is(g, BooleanGate):
            operands = list(dict.fromkeys(o.strip() for o in g.operands))
            if not operands:
                raise Unwritable("it is made of no population")
            if g.op == NOT and len(operands) > 1:
                any_ = GateNode(_unique(f"{name} [any]", taken), BooleanGate(OR, operands))
                any_.helper = True
                report.note(sample.name, node.path,
                            "the events in none of several populations are written as "
                            "the negation of their union, which is beside it")
                return [any_, GateNode(name, BooleanGate(NOT, [any_.name]))]
            if g.op in (AND, OR) and len(operands) == 1:
                operands = operands * 2      # both formats ask for two
            if g.op not in (AND, OR, NOT):
                raise Unwritable(f"a combination ({g.op}) that is not and, or, not")
            return [GateNode(name, BooleanGate(g.op, operands))]
        for kind in (RectangleGate, PolygonGate, EllipseGate, QuadrantGate):
            if _is(g, kind):
                return [GateNode(name, _checked(g))]
        raise Unwritable(f"a kind of gate ({type(g).__name__}) this export does not know")

    def run(src: GateNode, dst: GateNode) -> None:
        taken = {c.name.strip() for c in src.children}
        seen: set[str] = set()
        made: dict[int, list[GateNode]] = {}
        # the drawn gates first, as they are applied: a combination is
        # made of what is there by then
        drawn = [c for c in src.children if not isinstance(c.gate, BooleanGate)]
        for child in drawn + [c for c in src.children if isinstance(c.gate, BooleanGate)]:
            name = child.name.strip()
            try:
                if not name:
                    raise Unwritable("it has no name")
                if name in seen:
                    raise Unwritable("another population of the same parent has its name")
                if isinstance(child.gate, BooleanGate):
                    for o in child.gate.operands:
                        if o.strip() not in seen:
                            raise Unwritable(
                                f"it is made of {o.strip()}, which is not there to combine: "
                                "not beside it, not written, or a combination listed after it")
                made[id(child)] = expand(child, name, taken)
            except Unwritable as why:
                report.skip(sample.name, child, str(why))
                continue
            except Exception as exc:     # a gate this was not written for
                report.skip(sample.name, child, f"it could not be read: {exc}")
                continue
            seen.add(name)
        for child in src.children:
            if id(child) not in made:
                continue
            for node in made[id(child)]:
                dst.add(node)
            made[id(child)][-1].source = child.path
            run(child, made[id(child)][-1])

    run(sample.gates, root)
    return root


def counted(sample: ExportSample, tree: GateNode, report: Report) -> dict[str, int] | None:
    """The events of every population of the prepared tree, and a check on
    the way: a population written in another form must hold the events it
    holds here."""
    if sample.load is None:
        return None
    try:
        opened = sample.load()
        masks = tree.apply(opened.column, len(opened))
        own = opened.masks()
    except Exception as exc:
        report.note(sample.name, "", f"the counts are not written: {exc}")
        return None
    for node in tree.walk():
        path = getattr(node, "source", None)
        if node.gate is not None and path in own \
                and not np.array_equal(masks[node.path], own[path]):
            report.note(sample.name, path,
                        f"as written it holds {int(masks[node.path].sum()):,d} events, "
                        f"here {int(own[path].sum()):,d}")
    return {path: int(mask.sum()) for path, mask in masks.items()}


def walk(parent: GateNode, build, sample: str, report: Report) -> list[tuple[GateNode, object]]:
    """`build(node)` for every child of `parent` — the drawn gates first,
    so that a combination knows whether what it names was written — and
    (node, what it gave) for those that were, in the tree's order. A gate
    that cannot be built is left out and reported; so is a combination
    made of one."""
    written: set[str] = set()
    built: dict[int, object] = {}
    drawn = [c for c in parent.children if not isinstance(c.gate, BooleanGate)]
    for child in drawn + [c for c in parent.children if isinstance(c.gate, BooleanGate)]:
        try:
            if isinstance(child.gate, BooleanGate):
                for name in child.gate.operands:
                    if name not in written:
                        raise Unwritable(
                            f"it is made of {name}, which is not there to combine: not "
                            "beside it, not written, or a combination listed after it")
            built[id(child)] = build(child)
            written.add(child.name)
        except Unwritable as why:
            report.skip(sample, child, str(why))
        except Exception as exc:
            report.skip(sample, child, f"it could not be written: {exc}")
    return [(c, built[id(c)]) for c in parent.children if id(c) in built]


# --- scales and shapes -----------------------------------------------------------

def affine(scale) -> tuple[float, float] | None:
    """(a, b) where a stored value is a·s + b at the position s of the
    scale; None for a scale that bends."""
    if scale is None:
        return 1.0, 0.0
    if isinstance(scale, Linear):
        return scale.max - scale.min, scale.min
    return None


def on(scale, values):
    return np.asarray(values, dtype=np.float64) if scale is None else scale(values)


def off(scale, positions):
    return np.asarray(positions, dtype=np.float64) if scale is None else scale.inverse(positions)


def scales_of(gate) -> list:
    return list(getattr(gate, "scales", None) or [None, None])


def axis_scales(sample: ExportSample, tree: GateNode) -> dict[str, object]:
    """Gate dimension → the scale to give that parameter in the file. It
    is the scale of the axis in MIL-C — unless every polygon and ellipse
    drawn on the parameter was drawn on one other scale, which is then the
    one that keeps them as they are."""
    drawn: dict[str, list] = {}
    for node in tree.walk():
        if isinstance(node.gate, (PolygonGate, EllipseGate)):
            for dim, scale in zip(node.gate.dims, scales_of(node.gate)):
                drawn.setdefault(dim, []).append(scale)
    out = {}
    for dim in {d for n in tree.walk() if n.gate is not None for d in n.gate.dims}:
        usual = sample.scales.get(strip_prefix(dim)[0] if dim not in sample.channels else dim)
        used = drawn.get(dim, [])
        bent = [s for s in used if affine(s) is None]
        if bent and all(same(s, bent[0]) for s in used):
            out[dim] = bent[0]
        elif used and not bent and affine(usual) is None:
            out[dim] = None
        else:
            out[dim] = usual
    return out


def bend(points: np.ndarray, forward: list, back: list, tolerance: float = TOLERANCE,
         depth: int = 12) -> np.ndarray:
    """A polygon whose sides are straight where it was drawn, on other
    axes: `forward[k]` takes a coordinate there, `back[k]` returns it, None
    for an axis that stays as it is (or is only stretched). A side that
    would bend gets vertices added until the straight pieces are within
    `tolerance` of it — measured where it was drawn, against the polygon's
    size."""
    points = np.asarray(points, dtype=np.float64)

    def there(p):
        return np.column_stack([f(p[:, k]) if f is not None else p[:, k]
                                for k, f in enumerate(forward)])

    def here(p):
        return np.column_stack([f(p[:, k]) if f is not None else p[:, k]
                                for k, f in enumerate(back)])

    if all(f is None for f in forward):
        return points.copy()
    size = np.ptp(points, axis=0)
    size[size <= 0] = 1.0
    t = np.array([0.25, 0.5, 0.75])[None, :, None]
    for _ in range(depth):
        a, b = points, np.roll(points, -1, axis=0)
        fa = there(points)
        fb = np.roll(fa, -1, axis=0)
        # the straight piece between two vertices there, seen from here
        chord = here((fa[:, None, :] + (fb - fa)[:, None, :] * t).reshape(-1, 2)).reshape(-1, 3, 2)
        side = a[:, None, :] + (b - a)[:, None, :] * t
        split = np.flatnonzero((np.abs(chord - side) > tolerance * size).any(axis=(1, 2)))
        if not len(split):
            break
        points = np.insert(points, split + 1, (a[split] + b[split]) / 2, axis=0)
    return there(points)


def ellipse_points(center, radii, angle: float, sides: int = ELLIPSE_SIDES) -> np.ndarray:
    """An ellipse as a polygon of the same area, where it was drawn."""
    step = 2 * np.pi / sides
    grow = np.sqrt(step / np.sin(step))
    t = np.arange(sides) * step
    u, v = grow * radii[0] * np.cos(t), grow * radii[1] * np.sin(t)
    c, s = np.cos(angle), np.sin(angle)
    return np.column_stack([center[0] + u * c - v * s, center[1] + u * s + v * c])


def ellipse_covariance(radii, angle: float, stretch=(1.0, 1.0)) -> np.ndarray:
    """The matrix C of an ellipse (x−m)ᵀ C⁻¹ (x−m) ≤ 1, after each axis
    is multiplied by `stretch`."""
    c, s = np.cos(angle), np.sin(angle)
    turn = np.array([[c, -s], [s, c]])
    cov = turn @ np.diag([radii[0] ** 2, radii[1] ** 2]) @ turn.T
    a = np.diag(np.asarray(stretch, dtype=np.float64))
    return a @ cov @ a


def ellipse_axes(cov: np.ndarray) -> tuple[float, float, float]:
    """(long half-axis, short half-axis, angle of the long one)."""
    w, v = np.linalg.eigh(cov)
    return (float(np.sqrt(w[1])), float(np.sqrt(w[0])),
            float(np.arctan2(v[1, 1], v[0, 1])))


_fits: dict[tuple, tuple] = {}


def logicle_like(biex: Biex) -> tuple[Logicle, float] | None:
    """The Logicle scale nearest a FlowJo biexponential, and how far apart
    the two are at most between the ends of the biexponential, in the
    scale's own units (0 to 1; times 256 for FlowJo's channels). The two
    are the same curve but for how FlowJo builds its table, so they come
    within a few 1e-5. None where no Logicle is near."""
    key = (biex.length, biex.max_range, biex.neg, biex.width, biex.pos)
    if key not in _fits:
        _fits[key] = _fit_logicle(biex)
    return _fits[key]


def _fit_logicle(biex: Biex):
    from scipy.optimize import least_squares
    y = np.linspace(0.0, 1.0, 1025)
    x = biex.inverse(y)
    top = float(x[-1])
    if not np.all(np.diff(x) > 0) or top <= 0:
        return None
    zero = float(biex(0.0))
    # where the biexponential is nearly a Logicle: the decades of its
    # table, its zero at x1 = (W + A) / (M + A)
    width = np.log10(-biex.width) if biex.width < 0 else 0.0
    span = biex.pos + max(biex.neg, 0.0)
    start = np.array([span, width / 2, zero * span - width / 2])

    def make(p):
        m, w, a = p
        return Logicle(top, w, m - a, a)

    def apart(p):
        try:
            return make(p)(x) - y
        except (ValueError, FloatingPointError):
            return np.full(y.shape, 10.0)

    try:
        best = least_squares(apart, start, x_scale=[1.0, 0.1, 0.1], diff_step=1e-6)
        scale = make(best.x)
        worst = float(np.abs(scale(x) - y).max())
    except Exception:
        return None
    return (scale, worst) if worst < 0.01 else None


# --- from the session ----------------------------------------------------------------

def samples_of(session, entries=None, counts: bool = True) -> list[ExportSample]:
    """The session's samples as they are gated now: each with the matrix
    in use for it, the scales of the axes and its own tree."""
    from .fcs import read_header
    out = []
    for e in (entries if entries is not None else session.usable()):
        if e.problem:
            continue
        try:
            keywords = read_header(e.path)
        except Exception:
            keywords = {}
        matrix, source = session.matrix_in_use(e)
        out.append(ExportSample(
            name=os.path.basename(e.path), path=os.path.abspath(e.path), events=e.events,
            channels=list(e.channels), keywords=keywords, compensation=matrix,
            matrix_source=source, scales={c: session.scale_of(c) for c in e.channels},
            gates=e.gates, load=(lambda e=e: session.load(e)) if counts else None))
    return out


def export_gates(session, path: str, format: str | None = None, entries=None,
                 counts: bool = True) -> Report:
    """Write the gates of the session's samples (`entries`, or every usable
    one) to `path`, as a FlowJo workspace or as Gating-ML — by the file's
    extension where `format` is not given. What was left out is in the
    report that comes back."""
    from . import gatingml, wsp_write
    if format is None:
        format = WSP if path.lower().endswith(".wsp") else GATINGML
    samples = samples_of(session, entries, counts)
    report = (wsp_write.write_wsp if format == WSP else gatingml.write_gatingml)(samples, path)
    if session.normalization and session.normalization.get("everywhere") and report.written:
        report.note("", "", "the gates were drawn on batch-normalised values, which another "
                            "program does not have: it will count other events in them")
    return report


def spec_key(scale) -> str:
    return repr(to_spec(scale))
