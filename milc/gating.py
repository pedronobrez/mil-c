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
    gate: RectangleGate | PolygonGate | None = None
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

        def run(node: GateNode, parent_mask: np.ndarray):
            mask = parent_mask if node.gate is None else (
                parent_mask & node.gate.contains(column))
            out[node.path] = mask
            for c in node.children:
                run(c, mask)

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
        out = GateNode(node.name, g,
                       reference_count=node.reference_count if keep_counts else None)
        for child in node.children:
            out.add(clone(child))
        return out
    return clone(tree)


def _scale_out(scale):
    from .transforms import Biex
    if isinstance(scale, Biex):
        return {"biex": [scale.length, scale.max_range, scale.neg, scale.width, scale.pos]}
    return None


def _scale_in(spec):
    from .transforms import Biex
    return Biex(*spec["biex"]) if isinstance(spec, dict) and "biex" in spec else None


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
    node = GateNode(d.get("name", ""), gate)
    for child in d.get("children", []):
        node.add(tree_from_dict(child))
    return node
