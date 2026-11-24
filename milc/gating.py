"""
Gates and the tree they form.

A gate answers one question — which events are inside — on columns it asks
for by name. Where the columns come from (raw, compensated) is the
caller's business: `GateNode.apply` takes a `column(name)` function.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass, field
from typing import Callable, Iterator

import numpy as np

Column = Callable[[str], np.ndarray]


@dataclass
class RectangleGate:
    """One or two dimensions; a missing bound is open on that side."""
    dims: list[str]
    mins: list[float | None]
    maxs: list[float | None]
    inside: bool = True

    def contains(self, column: Column) -> np.ndarray:
        mask = None
        for d, lo, hi in zip(self.dims, self.mins, self.maxs):
            x = column(d)
            m = np.ones(x.shape[0], dtype=bool)
            if lo is not None:
                m &= x >= lo
            if hi is not None:
                m &= x < hi
            mask = m if mask is None else mask & m
        return mask if self.inside else ~mask


@dataclass
class PolygonGate:
    dims: list[str]
    #: n x 2, in the order of `dims`, in the units the file stores
    vertices: np.ndarray
    inside: bool = True
    #: per dimension, the scale the polygon was drawn on (a transform, or
    #: None for the values as they are). A side drawn straight on a
    #: biexponential axis is a curve in stored units, so the test is made
    #: where the polygon was drawn — which is also what FlowJo does.
    scales: list | None = None

    def contains(self, column: Column) -> np.ndarray:
        x, y = column(self.dims[0]), column(self.dims[1])
        poly = np.array(self.vertices, dtype=np.float64)
        for k, scale in enumerate(self.scales or ()):
            if scale is not None:
                poly[:, k] = scale(poly[:, k])
                if k == 0:
                    x = scale(x)
                else:
                    y = scale(y)
        mask = points_in_polygon(x, y, poly)
        return mask if self.inside else ~mask


def _scaled(values, scale):
    return scale(values) if scale is not None else values


@dataclass
class EllipseGate:
    """An ellipse on two parameters. Its centre, half-axes and rotation
    are kept in the coordinates it was drawn in — after `scales` — because
    an ellipse drawn on a biexponential axis is not an ellipse in the
    file's units, and it is the drawn one that is meant."""
    dims: list[str]
    center: list[float]
    radii: list[float]
    #: counter-clockwise, in radians, of the first half-axis from the X axis
    angle: float = 0.0
    inside: bool = True
    scales: list | None = None

    def contains(self, column: Column) -> np.ndarray:
        scales = self.scales or [None, None]
        dx = _scaled(column(self.dims[0]).astype(np.float64), scales[0]) - self.center[0]
        dy = _scaled(column(self.dims[1]).astype(np.float64), scales[1]) - self.center[1]
        c, s = np.cos(self.angle), np.sin(self.angle)
        u, v = dx * c + dy * s, -dx * s + dy * c
        mask = (u / self.radii[0]) ** 2 + (v / self.radii[1]) ** 2 <= 1.0
        return mask if self.inside else ~mask


#: the four corners a quadrant gate cuts the plot into, as FlowJo numbers
#: them: Q1 upper left, Q2 upper right, Q3 lower right, Q4 lower left.
#: Each is (X is high, Y is high).
QUADRANTS = {"Q1": (False, True), "Q2": (True, True), "Q3": (True, False), "Q4": (False, False)}


@dataclass
class QuadrantGate:
    """One corner of a plot cut in four at (x, y), in the file's units.
    The four corners made together share `group`, so moving the crossing
    point moves them all; each is still a population of its own."""
    dims: list[str]
    x: float
    y: float
    quadrant: str
    group: str = ""
    inside: bool = True

    def contains(self, column: Column) -> np.ndarray:
        high_x, high_y = QUADRANTS[self.quadrant]
        a, b = column(self.dims[0]) >= self.x, column(self.dims[1]) >= self.y
        mask = (a if high_x else ~a) & (b if high_y else ~b)
        return mask if self.inside else ~mask


#: the four arms of an offset quadrant, in the order they come going round
#: counter-clockwise, and the two that bound each corner
ARMS = ("right", "up", "left", "down")
SECTORS = {"Q2": ("right", "up"), "Q1": ("up", "left"),
           "Q4": ("left", "down"), "Q3": ("down", "right")}
_AXES = {"right": (1.0, 0.0), "up": (0.0, 1.0), "left": (-1.0, 0.0), "down": (0.0, -1.0)}


def _in_sector(a, b, dx: np.ndarray, dy: np.ndarray) -> np.ndarray:
    """Is (dx, dy) in the closed sector swept counter-clockwise from the
    direction `a` to the direction `b`?"""
    after_a = a[0] * dy - a[1] * dx >= 0
    before_b = dx * b[1] - dy * b[0] >= 0
    turn = a[0] * b[1] - a[1] * b[0]
    if turn > 0:                        # less than half a turn
        return after_a & before_b
    if turn < 0:                        # more than half a turn
        return after_a | before_b
    if a[0] * b[0] + a[1] * b[1] < 0:   # exactly half: one side of a line
        return after_a
    # the two arms on top of each other: nothing between them but the arm
    return after_a & before_b & (a[0] * dx + a[1] * dy >= 0)


@dataclass
class OffsetQuadrantGate:
    """One corner of a plot cut in four by arms that need not be at right
    angles nor in line — FlowJo's quadrant with its arms offset.

    The geometry: a `center` and four `arms`, `right`, `up`, `left` and
    `down`, each given by one point it passes through. An arm is the
    half-line that starts at the centre, goes through its point and on
    for ever. Going round counter-clockwise the arms come in that order,
    and each corner is the sector between two neighbours: Q2 (upper
    right) from `right` to `up`, Q1 (upper left) from `up` to `left`, Q4
    (lower left) from `left` to `down`, Q3 (lower right) from `down` to
    `right`. An event on an arm goes where `QuadrantGate` puts it: on
    `right` or `up` in Q2, on `left` in Q1, on `down` in Q3, and at the
    centre in Q2 — so with the arms along the axes this is `QuadrantGate`
    at the centre, event for event. Every event is in exactly one corner,
    whatever the arms.

    `center` and the arms' points are in the file's units, in the order
    of `dims`; the arms are straight on `scales` — the display scale of
    each axis, None for the values as they are — which is where the
    events are tested, as for `PolygonGate`. An arm whose point is the
    centre has no direction and is taken along its axis. The four corners
    made together share `group`, the centre and the arms."""
    dims: list[str]
    center: list[float]
    #: "right", "up", "left", "down" → [x, y], a point the arm goes through
    arms: dict[str, list[float]]
    quadrant: str
    group: str = ""
    inside: bool = True
    scales: list | None = None

    def shown(self) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        """The centre and each arm's point on the display scales."""
        scales = self.scales or [None, None]

        def place(p):
            return np.array([float(_scaled(np.float64(v), s)) for v, s in zip(p, scales)])
        return place(self.center), {arm: place(self.arms[arm]) for arm in ARMS}

    def directions(self) -> dict[str, tuple[float, float]]:
        """Each arm's direction from the centre, on the display scales."""
        centre, points = self.shown()
        out = {}
        for arm in ARMS:
            d = points[arm] - centre
            out[arm] = ((float(d[0]), float(d[1]))
                        if np.all(np.isfinite(d)) and np.any(d != 0) else _AXES[arm])
        return out

    def ordered(self) -> bool:
        """Do the arms come in their order going round, none on or past
        its neighbour? If not, the corners still share every event
        between them, but are not the four their names say."""
        d = self.directions()
        angles = [np.arctan2(d[a][1], d[a][0]) for a in ARMS]
        steps = [(angles[(k + 1) % 4] - angles[k]) % (2 * np.pi) for k in range(4)]
        return all(s > 1e-9 for s in steps) and abs(sum(steps) - 2 * np.pi) < 1e-6

    def corner(self, xr: tuple[float, float], yr: tuple[float, float]) -> np.ndarray | None:
        """This corner as a polygon, n x 2 on the display scales: its
        sector cut off at the box from `xr` to `yr` — for a program that
        has polygons and no such quadrant. The centre first, then round
        counter-clockwise. None when the centre is not inside the box."""
        centre, _ = self.shown()
        if not (xr[0] < centre[0] < xr[1] and yr[0] < centre[1] < yr[1]):
            return None
        d = self.directions()
        first, second = (np.array(d[arm]) for arm in SECTORS[self.quadrant])

        def leaves(v):
            steps = [((hi if u > 0 else lo) - c) / u
                     for u, c, (lo, hi) in zip(v, centre, (xr, yr)) if u != 0]
            return centre + min(steps) * v

        def angle(p):
            return np.arctan2(p[1] - centre[1], p[0] - centre[0])
        start = angle(centre + first)
        sweep = (angle(centre + second) - start) % (2 * np.pi)
        corners = [np.array(p, dtype=np.float64)
                   for p in ((xr[1], yr[1]), (xr[0], yr[1]), (xr[0], yr[0]), (xr[1], yr[0]))]
        between = sorted((((angle(p) - start) % (2 * np.pi), k)
                          for k, p in enumerate(corners)))
        return np.array([centre, leaves(first)]
                        + [corners[k] for turn, k in between if 0 < turn < sweep]
                        + [leaves(second)])

    def contains(self, column: Column) -> np.ndarray:
        scales = self.scales or [None, None]
        centre, _ = self.shown()
        dx = _scaled(column(self.dims[0]).astype(np.float64), scales[0]) - centre[0]
        dy = _scaled(column(self.dims[1]).astype(np.float64), scales[1]) - centre[1]
        d = self.directions()
        taken = np.zeros(dx.shape[0], dtype=bool)
        mask = taken
        # in this order an event on an arm, which two sectors hold, goes
        # to the corner it belongs to; Q4 is what the other three leave
        for quadrant in ("Q2", "Q1", "Q3", "Q4"):
            a, b = SECTORS[quadrant]
            mask = ~taken if quadrant == "Q4" else _in_sector(d[a], d[b], dx, dy) & ~taken
            if quadrant == self.quadrant:
                break
            taken = taken | mask
        return mask if self.inside else ~mask


def offset_quadrant(gate: "QuadrantGate", scales: list | None = None,
                    reach: tuple[float, float] = (1.0, 1.0)) -> OffsetQuadrantGate:
    """The same corner with arms that can be moved: along the axes, so
    not one event changes side until an arm is moved. `reach` is how far
    from the centre, on the display scales, each arm's point is put."""
    scales = list(scales) if scales and any(s is not None for s in scales) else None
    sc = scales or [None, None]
    cx, cy = (float(_scaled(np.float64(v), s)) for v, s in zip((gate.x, gate.y), sc))

    def stored(value, k):
        return float(sc[k].inverse(np.float64(value))) if sc[k] is not None else float(value)
    arms = {}
    for arm, (ux, uy) in _AXES.items():
        # the coordinate an arm does not move along stays the centre's own
        # number, not one taken there and back through the scale
        arms[arm] = [stored(cx + ux * reach[0], 0) if ux else float(gate.x),
                     stored(cy + uy * reach[1], 1) if uy else float(gate.y)]
    return OffsetQuadrantGate(list(gate.dims), [float(gate.x), float(gate.y)], arms,
                              gate.quadrant, gate.group, gate.inside, scales)


#: the curve of a curly quadrant as CytoML draws it (cytolib, src/gate.cpp,
#: CurlyQuadGate::interpolate): on 255 channels of each axis, the arm to
#: the right rises by 0.001·d², d the channels from the centre, and the arm
#: going up leans right by as much. FlowJo writes no curve in a workspace.
CURL_CHANNELS = 255.0
CURL_POWER = 2.0
CURL_BEND = 0.001
#: where cytolib puts the top of an axis FlowJo gives no scale of its own
CURL_TOP = 262144.0


def _unbounded(scale) -> bool:
    """Does the scale go on past its ends, so that a point beyond them is
    still told apart? FlowJo's biexponential is a table, and stops."""
    from .transforms import Biex
    return not isinstance(scale, Biex)


@dataclass
class CurlyQuadrantGate:
    """One corner of a plot cut in four by FlowJo's "curly quadrant": the
    arms to the left of and below the centre are straight, the two others
    bend away from it, the way the spread of a compensated negative
    population grows with the intensity on the other axis.

    The curve is measured in channels: each axis' position on its scale
    (`scales`, from 0 to 1) times `channels`; an axis with no scale is
    linear from 0 to `CURL_TOP`. With (cx, cy) the centre there and
    (x, y) an event, it is

        above the horizontal divider:  y >= cy + bend · max(x - cx, 0) ** power
        right of the vertical divider: x >= cx + bend · max(y - cy, 0) ** power

    and Q2 (upper right) is both, Q1 (upper left) only the first, Q3 (lower
    right) only the second, Q4 neither — so every event is in exactly one
    corner, one on a divider going to the side above or to the right, as
    for `QuadrantGate`, which this is with `bend` 0.

    `center` is in the file's units, in the order of `dims`. The default
    curve is CytoML's (`CURL_BEND`, `CURL_POWER` on `CURL_CHANNELS`): FlowJo
    computes its own and writes none of it (see milc/wsp.py)."""
    dims: list[str]
    center: list[float]
    quadrant: str
    group: str = ""
    inside: bool = True
    scales: list | None = None
    bend: float = CURL_BEND
    power: float = CURL_POWER
    channels: float = CURL_CHANNELS

    # --- where the curve is measured -------------------------------------------

    def _per(self, k: int) -> float:
        """Channels in one unit of the axis' shown coordinate."""
        scale = (self.scales or [None, None])[k]
        return self.channels if scale is not None else self.channels / CURL_TOP

    def to_channels(self, values, k: int) -> np.ndarray:
        scale = (self.scales or [None, None])[k]
        return _scaled(np.asarray(values, dtype=np.float64), scale) * self._per(k)

    def from_channels(self, channels, k: int) -> np.ndarray:
        scale = (self.scales or [None, None])[k]
        shown = np.asarray(channels, dtype=np.float64) / self._per(k)
        return scale.inverse(shown) if scale is not None else shown

    def centre_channels(self) -> tuple[float, float]:
        return tuple(float(self.to_channels(np.float64(v), k))
                     for k, v in enumerate(self.center))

    def rise(self, d) -> np.ndarray:
        """How far an arm has bent `d` channels from the centre."""
        d = np.maximum(np.asarray(d, dtype=np.float64), 0.0)
        return self.bend * d ** self.power if self.bend else np.zeros_like(d)

    def sides(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(right of the vertical divider, above the horizontal one) for
        events at (x, y) in channels."""
        cx, cy = self.centre_channels()
        with np.errstate(invalid="ignore", over="ignore"):
            right = x >= cx + self.rise(y - cy)
            up = y >= cy + self.rise(x - cx)
        return right, up

    def contains(self, column: Column) -> np.ndarray:
        x = self.to_channels(column(self.dims[0]), 0)
        y = self.to_channels(column(self.dims[1]), 1)
        right, up = self.sides(x, y)
        high_x, high_y = QUADRANTS[self.quadrant]
        mask = (right if high_x else ~right) & (up if high_y else ~up)
        return mask if self.inside else ~mask

    # --- its shape ----------------------------------------------------------------

    def arm_channels(self, arm: str, reach: float, steps: int = 96) -> np.ndarray:
        """One arm from the centre out to `reach` channels along its axis,
        n x 2 in channels: bent for `right` and `up`, straight for the
        others."""
        cx, cy = self.centre_channels()
        d = np.linspace(0.0, max(float(reach), 0.0), steps if arm in ("right", "up") else 2)
        if arm == "right":
            return np.column_stack([cx + d, cy + self.rise(d)])
        if arm == "up":
            return np.column_stack([cx + self.rise(d), cy + d])
        if arm == "left":
            return np.column_stack([cx - d, np.full_like(d, cy)])
        return np.column_stack([np.full_like(d, cx), cy - d])

    def arms(self, reach: float | None = None, steps: int = 96) -> dict[str, np.ndarray]:
        """The four arms in the file's units, each from the centre to
        `reach` channels along its axis (by default twice the scale) — a
        bent one only as far as both its coordinates are on their scales:
        past the end of a scale that stops there (a biexponential) a point
        cannot be told from the end, and the curve would run along it."""
        reach = 2 * self.channels if reach is None else reach
        centre = self.centre_channels()
        tops = [self.channels if s is not None and not _unbounded(s) else np.inf
                for s in (self.scales or [None, None])]
        out = {}
        for arm in ARMS:
            far = reach
            if arm in ("right", "up") and self.bend > 0:
                along, across = (0, 1) if arm == "right" else (1, 0)
                far = min(far, tops[along] - centre[along])
                room = tops[across] - centre[across]
                if np.isfinite(room) and room > 0:
                    far = min(far, (room / self.bend) ** (1.0 / self.power))
            pts = self.arm_channels(arm, max(far, 0.0), steps)
            out[arm] = np.column_stack([self.from_channels(pts[:, k], k) for k in range(2)])
        return out

    def corner(self, xr: tuple[float, float], yr: tuple[float, float],
               steps: int = 256) -> np.ndarray | None:
        """This corner as a polygon, n x 2 on the gate's scales (the file's
        units on an axis with none): bounded by its two arms, each followed
        in `steps` straight pieces, and the box from `xr` to `yr`. The
        centre first, then round counter-clockwise. None when the centre is
        not inside the box."""
        per = [self._per(0), self._per(1)]
        box = [(lo * p, hi * p) for (lo, hi), p in zip((xr, yr), per)]
        cx, cy = self.centre_channels()
        (x0, x1), (y0, y1) = box
        if not (x0 < cx < x1 and y0 < cy < y1):
            return None

        def inside(p):
            return (p[:, 0] >= x0) & (p[:, 0] <= x1) & (p[:, 1] >= y0) & (p[:, 1] <= y1)

        def clipped(arm):
            """The arm up to where it leaves the box, that point included."""
            reach = max(x1 - x0, y1 - y0)
            pts = self.arm_channels(arm, reach, steps + 1)
            while bool(inside(pts[-1:])[0]):     # not out yet: further
                reach *= 2
                pts = self.arm_channels(arm, reach, steps + 1)
            out_at = int(np.argmin(inside(pts)))
            a, b = pts[out_at - 1], pts[out_at]
            # where the piece from the last point in to the first out
            # crosses the box
            t = 1.0
            for k, (lo, hi) in enumerate(box):
                for edge in (lo, hi):
                    if (b[k] - edge) * (a[k] - edge) < 0:
                        t = min(t, (edge - a[k]) / (b[k] - a[k]))
            return np.vstack([pts[:out_at], a + t * (b - a)])

        first, second = (clipped(arm) for arm in SECTORS[self.quadrant])
        centre = np.array([cx, cy])

        def angle(p):
            return np.arctan2(p[1] - cy, p[0] - cx)
        start = angle(first[-1])
        sweep = (angle(second[-1]) - start) % (2 * np.pi)
        corners = [np.array(p, dtype=np.float64) for p in ((x1, y1), (x0, y1), (x0, y0), (x1, y0))]
        between = sorted((((angle(p) - start) % (2 * np.pi), k) for k, p in enumerate(corners)))
        # the box is seen from the centre one turn round, so its corners
        # between where the two arms leave it are those between their angles
        ring = np.vstack([centre, first[1:]]
                         + [corners[k][None, :] for turn, k in between if 0 < turn < sweep]
                         + [second[:0:-1]])
        return ring / np.array(per)


AND, OR, NOT = "and", "or", "not"


@dataclass
class BooleanGate:
    """A population made of its sisters: the events in all of them, in
    any of them, or in none of them. `operands` are the names of other
    populations of the same parent. It has no place on a plot."""
    op: str
    operands: list[str]
    dims: list[str] = field(default_factory=list)
    inside: bool = True

    def combine(self, masks: dict[str, np.ndarray], n: int) -> np.ndarray:
        found = [masks[name] for name in self.operands if name in masks]
        if len(found) < len(self.operands) or not found:
            # something it is made of is gone: nothing, rather than a guess
            return np.zeros(n, dtype=bool)
        if self.op == AND:
            return np.logical_and.reduce(found)
        if self.op == OR:
            return np.logical_or.reduce(found)
        return ~np.logical_or.reduce(found)


@dataclass
class ClusterGate:
    """The events a trained self-organising map puts in one cluster: each
    event goes to its nearest node, on the channels and scales the map was
    trained on, and is kept if that node is one of the cluster's. It
    carries the map, so it means the same on any file and in any project
    it is saved in. It has no place on a plot."""
    #: the columns the map was trained on, `Comp-` where compensated
    channels: list[str]
    #: per channel, the scale it was given (None: as stored) and what that
    #: was multiplied by on the way to the map
    scales: list
    factors: list[float]
    #: nodes x channels: the map
    codes: np.ndarray
    #: the nodes that belong to the cluster
    nodes: list[int]
    dims: list[str] = field(default_factory=list)
    inside: bool = True

    def contains(self, column: Column) -> np.ndarray:
        from .flowsom import _nearest
        # through float32, as the concatenated events are kept: an event
        # on the border between two nodes must fall here where it fell there
        X = np.column_stack([
            _scaled(column(c).astype(np.float32).astype(np.float64), s) * f
            for c, s, f in zip(self.channels, self.scales, self.factors)])
        mask = np.isin(_nearest(np.ascontiguousarray(X), self.codes), self.nodes)
        return mask if self.inside else ~mask


@dataclass
class TimeGate:
    """The events acquired inside — or, with `inside` off, outside — a list
    of time intervals, both ends included, in the units the file stores
    Time in. It is how quality control removes the stretches of an
    acquisition that went wrong, and it belongs to one sample: another
    sample's bad stretches are its own (see `keep_times`)."""
    dims: list[str]
    intervals: list[list[float]]
    inside: bool = True

    def contains(self, column: Column) -> np.ndarray:
        t = column(self.dims[0])
        mask = np.zeros(t.shape[0], dtype=bool)
        for lo, hi in self.intervals:
            mask |= (t >= lo) & (t <= hi)
        return mask if self.inside else ~mask


def points_in_polygon(x: np.ndarray, y: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Even-odd rule, one pass per edge."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    inside = np.zeros(x.shape[0], dtype=bool)
    # nothing outside the bounding box needs the edges
    box = ((x >= poly[:, 0].min()) & (x <= poly[:, 0].max())
           & (y >= poly[:, 1].min()) & (y <= poly[:, 1].max()))
    xb, yb = x[box], y[box]
    hit = np.zeros(xb.shape[0], dtype=bool)
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if y1 == y2:
            continue
        crosses = (y1 > yb) != (y2 > yb)
        xi = x1 + (yb - y1) * (x2 - x1) / (y2 - y1)
        hit ^= crosses & (xb < xi)
    inside[box] = hit
    return inside


@dataclass
class GateNode:
    """A population: a gate applied to its parent's events. The root has no
    gate and holds every event."""
    name: str
    gate: (RectangleGate | PolygonGate | EllipseGate | QuadrantGate
           | OffsetQuadrantGate | CurlyQuadrantGate | BooleanGate | None) = None
    children: list["GateNode"] = field(default_factory=list)
    #: the count another program reported for this population, when the
    #: tree was imported: what ours is checked against
    reference_count: int | None = None
    parent: "GateNode | None" = field(default=None, repr=False, compare=False)
    #: what a user should be told about where this gate came from, when
    #: it was imported by a rule not checked against the program that
    #: wrote it; empty when there is nothing to say
    note: str = ""
    #: how the gate was fitted to this sample's own events, when each
    #: sample was given its own (`milc.autogate`): the method and its
    #: settings, as plain values. A copy of the tree given to another
    #: sample then keeps that sample's own gate (see `keep_fits`). With a
    #: `sample` key the gate is not each sample's own: it was fitted on
    #: that one file's events and every sample is given the same.
    fit: dict | None = None
    #: what the fitted gate was fitted on, and what it came out as: the
    #: digests `on` (of the parent's gates and the values under them) and
    #: `gate` — see `Session.refit`. This sample's own, like the gate.
    fitted: dict | None = None

    def add(self, child: "GateNode") -> "GateNode":
        child.parent = self
        self.children.append(child)
        return child

    @property
    def path(self) -> str:
        """`Cells/Single Cells/Live`. Names repeat (two `Single Cells` in a
        row is ordinary), so a path, not a name, identifies a population."""
        parts, node = [], self
        while node is not None and node.gate is not None:
            parts.append(node.name)
            node = node.parent
        return "/".join(reversed(parts))

    def walk(self) -> Iterator["GateNode"]:
        yield self
        for c in self.children:
            yield from c.walk()

    def find(self, path: str) -> "GateNode":
        for node in self.walk():
            if node.path == path:
                return node
        raise KeyError(path)

    def leaves(self) -> list["GateNode"]:
        return [n for n in self.walk() if not n.children and n.gate is not None]

    def apply(self, column: Column, n_events: int) -> dict[str, np.ndarray]:
        """path → boolean mask over all the sample's events, for every
        population in the tree. The root is under ``""``."""
        out: dict[str, np.ndarray] = {}

        def run(node: GateNode, mask: np.ndarray):
            out[node.path] = mask
            sisters: dict[str, np.ndarray] = {}
            # the drawn gates first: a combination is made of them, in
            # whatever order the tree lists them
            ordered = ([c for c in node.children if not isinstance(c.gate, BooleanGate)]
                       + [c for c in node.children if isinstance(c.gate, BooleanGate)])
            for c in ordered:
                if isinstance(c.gate, BooleanGate):
                    inner = mask & c.gate.combine(sisters, n_events)
                else:
                    inner = mask & c.gate.contains(column)
                sisters[c.name] = inner
                run(c, inner)

        run(self, np.ones(n_events, dtype=bool))
        return out


# --- copies and files --------------------------------------------------------

def copy_tree(tree: GateNode, keep_counts: bool = True) -> GateNode:
    """An independent copy: editing it leaves the original alone."""
    def clone(node: GateNode) -> GateNode:
        g = node.gate
        if isinstance(g, RectangleGate):
            g = RectangleGate(list(g.dims), list(g.mins), list(g.maxs), g.inside)
        elif isinstance(g, PolygonGate):
            g = PolygonGate(list(g.dims), np.array(g.vertices, dtype=np.float64),
                            g.inside, list(g.scales) if g.scales else None)
        elif isinstance(g, EllipseGate):
            g = EllipseGate(list(g.dims), list(g.center), list(g.radii), g.angle,
                            g.inside, list(g.scales) if g.scales else None)
        elif isinstance(g, QuadrantGate):
            g = QuadrantGate(list(g.dims), g.x, g.y, g.quadrant, g.group, g.inside)
        elif isinstance(g, OffsetQuadrantGate):
            g = OffsetQuadrantGate(list(g.dims), list(g.center),
                                   {a: list(g.arms[a]) for a in ARMS}, g.quadrant,
                                   g.group, g.inside, list(g.scales) if g.scales else None)
        elif isinstance(g, CurlyQuadrantGate):
            g = CurlyQuadrantGate(list(g.dims), list(g.center), g.quadrant, g.group, g.inside,
                                  list(g.scales) if g.scales else None, g.bend, g.power,
                                  g.channels)
        elif isinstance(g, BooleanGate):
            g = BooleanGate(g.op, list(g.operands), inside=g.inside)
        elif isinstance(g, TimeGate):
            g = TimeGate(list(g.dims), [list(i) for i in g.intervals], g.inside)
        elif isinstance(g, ClusterGate):
            g = ClusterGate(list(g.channels), list(g.scales), list(g.factors),
                            np.array(g.codes, dtype=np.float64), list(g.nodes),
                            inside=g.inside)
        out = GateNode(node.name, g, note=node.note,
                       reference_count=node.reference_count if keep_counts else None,
                       fit=_plain_copy(node.fit), fitted=_plain_copy(node.fitted))
        for child in node.children:
            out.add(clone(child))
        return out
    return clone(tree)


def _plain_copy(value):
    if isinstance(value, dict):
        return {k: _plain_copy(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain_copy(v) for v in value]
    return value


def _scale_out(scale):
    from .transforms import to_spec
    return to_spec(scale)


def _scale_in(spec):
    from .transforms import from_spec
    return from_spec(spec)


def tree_to_dict(tree: GateNode) -> dict:
    out: dict = {"name": tree.name, "children": [tree_to_dict(c) for c in tree.children]}
    if tree.note:
        out["note"] = tree.note
    if tree.fit:
        out["fit"] = _plain_copy(tree.fit)
    if tree.fitted:
        out["fitted"] = _plain_copy(tree.fitted)
    g = tree.gate
    if isinstance(g, RectangleGate):
        out["gate"] = {"type": "rectangle", "dims": g.dims, "mins": g.mins,
                       "maxs": g.maxs, "inside": g.inside}
    elif isinstance(g, PolygonGate):
        out["gate"] = {"type": "polygon", "dims": g.dims, "inside": g.inside,
                       "vertices": np.asarray(g.vertices, dtype=float).tolist(),
                       "scales": [_scale_out(s) for s in (g.scales or [None, None])]}
    elif isinstance(g, EllipseGate):
        out["gate"] = {"type": "ellipse", "dims": g.dims, "inside": g.inside,
                       "center": [float(v) for v in g.center],
                       "radii": [float(v) for v in g.radii], "angle": float(g.angle),
                       "scales": [_scale_out(s) for s in (g.scales or [None, None])]}
    elif isinstance(g, QuadrantGate):
        out["gate"] = {"type": "quadrant", "dims": g.dims, "inside": g.inside,
                       "x": float(g.x), "y": float(g.y), "quadrant": g.quadrant,
                       "group": g.group}
    elif isinstance(g, OffsetQuadrantGate):
        out["gate"] = {"type": "offset_quadrant", "dims": g.dims, "inside": g.inside,
                       "center": [float(v) for v in g.center],
                       "arms": {a: [float(v) for v in g.arms[a]] for a in ARMS},
                       "quadrant": g.quadrant, "group": g.group,
                       "scales": [_scale_out(s) for s in (g.scales or [None, None])]}
    elif isinstance(g, CurlyQuadrantGate):
        out["gate"] = {"type": "curly_quadrant", "dims": g.dims, "inside": g.inside,
                       "center": [float(v) for v in g.center], "quadrant": g.quadrant,
                       "group": g.group, "bend": float(g.bend), "power": float(g.power),
                       "channels": float(g.channels),
                       "scales": [_scale_out(s) for s in (g.scales or [None, None])]}
    elif isinstance(g, BooleanGate):
        out["gate"] = {"type": "boolean", "op": g.op, "operands": list(g.operands)}
    elif isinstance(g, ClusterGate):
        # the map as its bytes: a tenth of the size of its digits, and exact
        codes = np.ascontiguousarray(g.codes, dtype="<f8")
        out["gate"] = {"type": "cluster", "channels": list(g.channels),
                       "scales": [_scale_out(s) for s in g.scales],
                       "factors": [float(f) for f in g.factors],
                       "nodes": [int(n) for n in g.nodes], "inside": g.inside,
                       "codes": {"shape": list(codes.shape),
                                 "float64": base64.b64encode(codes.tobytes()).decode("ascii")}}
    elif isinstance(g, TimeGate):
        out["gate"] = {"type": "time", "dims": g.dims, "inside": g.inside,
                       "intervals": [[float(lo), float(hi)] for lo, hi in g.intervals]}
    return out


def tree_from_dict(d: dict) -> GateNode:
    g = d.get("gate")
    gate = None
    if g and g["type"] == "rectangle":
        gate = RectangleGate(list(g["dims"]), list(g["mins"]), list(g["maxs"]),
                             g.get("inside", True))
    elif g and g["type"] == "polygon":
        scales = [_scale_in(s) for s in g.get("scales") or []]
        gate = PolygonGate(list(g["dims"]), np.array(g["vertices"], dtype=np.float64),
                           g.get("inside", True),
                           scales if any(s is not None for s in scales) else None)
    elif g and g["type"] == "ellipse":
        scales = [_scale_in(s) for s in g.get("scales") or []]
        gate = EllipseGate(list(g["dims"]), list(g["center"]), list(g["radii"]),
                           float(g.get("angle", 0.0)), g.get("inside", True),
                           scales if any(s is not None for s in scales) else None)
    elif g and g["type"] == "quadrant":
        gate = QuadrantGate(list(g["dims"]), float(g["x"]), float(g["y"]), g["quadrant"],
                            g.get("group", ""), g.get("inside", True))
    elif g and g["type"] == "offset_quadrant":
        scales = [_scale_in(s) for s in g.get("scales") or []]
        gate = OffsetQuadrantGate(list(g["dims"]), [float(v) for v in g["center"]],
                                  {a: [float(v) for v in g["arms"][a]] for a in ARMS},
                                  g["quadrant"], g.get("group", ""), g.get("inside", True),
                                  scales if any(s is not None for s in scales) else None)
    elif g and g["type"] == "curly_quadrant":
        scales = [_scale_in(s) for s in g.get("scales") or []]
        gate = CurlyQuadrantGate(list(g["dims"]), [float(v) for v in g["center"]],
                                 g["quadrant"], g.get("group", ""), g.get("inside", True),
                                 scales if any(s is not None for s in scales) else None,
                                 float(g.get("bend", CURL_BEND)), float(g.get("power", CURL_POWER)),
                                 float(g.get("channels", CURL_CHANNELS)))
    elif g and g["type"] == "boolean":
        gate = BooleanGate(g["op"], list(g["operands"]))
    elif g and g["type"] == "time":
        gate = TimeGate(list(g["dims"]), [[float(lo), float(hi)] for lo, hi in g["intervals"]],
                        g.get("inside", True))
    elif g and g["type"] == "cluster":
        codes = np.frombuffer(base64.b64decode(g["codes"]["float64"]), dtype="<f8")
        gate = ClusterGate(list(g["channels"]), [_scale_in(s) for s in g["scales"]],
                           [float(f) for f in g["factors"]],
                           codes.reshape(g["codes"]["shape"]).astype(np.float64),
                           [int(n) for n in g["nodes"]], inside=g.get("inside", True))
    node = GateNode(d.get("name", ""), gate, note=d.get("note", ""),
                    fit=d.get("fit") or None, fitted=d.get("fitted") or None)
    for child in d.get("children", []):
        node.add(tree_from_dict(child))
    return node


def depends_on(node: GateNode, name: str) -> bool:
    """Is `node` a combination that names the sister `name`?"""
    return isinstance(node.gate, BooleanGate) and name in node.gate.operands


def follow_scale(tree: GateNode, dim: str, scale) -> bool:
    """A parameter's axis was given another scale: the polygons drawn on
    it are tested on the new one, as their sides are now straight there,
    and so are the arms of an offset quadrant, and a curly quadrant bends
    in channels of the new one, as FlowJo redraws it on the axes it is on.
    (An ellipse keeps the scale it was drawn on; it cannot be carried
    over without becoming another shape.) Says whether anything changed."""
    changed = False
    for node in tree.walk():
        g = node.gate
        if isinstance(g, (PolygonGate, OffsetQuadrantGate, CurlyQuadrantGate)) and dim in g.dims:
            scales = list(g.scales or [None, None])
            scales[g.dims.index(dim)] = scale
            g.scales = scales if any(x is not None for x in scales) else None
            changed = True
    return changed


# --- a gate put above every other: the acquisition's quality -----------------

def time_node(tree: GateNode) -> GateNode | None:
    """The population quality control made at the top of this tree."""
    return next((c for c in tree.children if isinstance(c.gate, TimeGate)), None)


def wrap(tree: GateNode, name: str, gate) -> GateNode:
    """Put a population at the top of the tree and everything that was
    there under it, so every later gate starts from its events."""
    node = GateNode(name, gate)
    below, tree.children = tree.children, []
    tree.add(node)
    for child in below:
        node.add(child)
    return node


def unwrap(tree: GateNode, node: GateNode) -> None:
    """Take a population of the root out again; what was under it goes
    back to where it was."""
    at = tree.children.index(node)
    tree.children[at:at + 1] = node.children
    for child in node.children:
        child.parent = tree
    node.children, node.parent = [], None


def keep_times(copy: GateNode, own: GateNode) -> None:
    """`copy` is another sample's tree, given to this one. Its time gates
    take the intervals this sample's `own` tree had — matched in order,
    as a gate may have been renamed — or none: the stretches of one
    acquisition that went wrong say nothing about another's."""
    had = [n.gate.intervals for n in own.walk() if isinstance(n.gate, TimeGate)]
    new = [n.gate for n in copy.walk() if isinstance(n.gate, TimeGate)]
    for k, gate in enumerate(new):
        gate.intervals = [list(i) for i in had[k]] if k < len(had) else []


def keep_fits(copy: GateNode, own: GateNode) -> list[GateNode]:
    """`copy` is another sample's tree, given to this one. A gate fitted
    on each sample's own events takes back the one this sample's `own`
    tree had at the same path, if that was fitted the same way. The
    populations it had none for are returned: their gates are still the
    other sample's, to be fitted here. A gate fitted on one sample and
    given to all (`fit` with a `sample`) is the copy's, as a drawn one is."""
    had = {n.path: n for n in own.walk() if n.fit}
    missing = []
    for node in copy.walk():
        if not node.fit or node.fit.get("sample"):
            continue
        mine = had.get(node.path)
        if mine is not None and mine.fit == node.fit:
            node.gate, node.note = copy_tree(GateNode("", mine.gate)).gate, mine.note
            node.fitted = _plain_copy(mine.fitted)
        else:
            missing.append(node)
    return missing
