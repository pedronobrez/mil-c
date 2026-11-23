"""
Gates and the tree they form.

A gate answers one question — which events are inside — on columns it asks
for by name. Where the columns come from (raw, compensated) is the
caller's business: `GateNode.apply` takes a `column(name)` function.
"""

from __future__ import annotations

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
           | BooleanGate | None) = None
    children: list["GateNode"] = field(default_factory=list)
    #: the count another program reported for this population, when the
    #: tree was imported: what ours is checked against
    reference_count: int | None = None
    parent: "GateNode | None" = field(default=None, repr=False, compare=False)

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
        elif isinstance(g, BooleanGate):
            g = BooleanGate(g.op, list(g.operands), inside=g.inside)
        out = GateNode(node.name, g,
                       reference_count=node.reference_count if keep_counts else None)
        for child in node.children:
            out.add(clone(child))
        return out
    return clone(tree)


def _scale_out(scale):
    from .transforms import to_spec
    return to_spec(scale)


def _scale_in(spec):
    from .transforms import from_spec
    return from_spec(spec)


def tree_to_dict(tree: GateNode) -> dict:
    out: dict = {"name": tree.name, "children": [tree_to_dict(c) for c in tree.children]}
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
    elif isinstance(g, BooleanGate):
        out["gate"] = {"type": "boolean", "op": g.op, "operands": list(g.operands)}
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
    elif g and g["type"] == "boolean":
        gate = BooleanGate(g["op"], list(g["operands"]))
    node = GateNode(d.get("name", ""), gate)
    for child in d.get("children", []):
        node.add(tree_from_dict(child))
    return node


def depends_on(node: GateNode, name: str) -> bool:
    """Is `node` a combination that names the sister `name`?"""
    return isinstance(node.gate, BooleanGate) and name in node.gate.operands


def follow_scale(tree: GateNode, dim: str, scale) -> bool:
    """A parameter's axis was given another scale: the polygons drawn on
    it are tested on the new one, as their sides are now straight there.
    (An ellipse keeps the scale it was drawn on; it cannot be carried
    over without becoming another shape.) Says whether anything changed."""
    changed = False
    for node in tree.walk():
        g = node.gate
        if isinstance(g, PolygonGate) and dim in g.dims:
            scales = list(g.scales or [None, None])
            scales[g.dims.index(dim)] = scale
            g.scales = scales if any(x is not None for x in scales) else None
            changed = True
    return changed
