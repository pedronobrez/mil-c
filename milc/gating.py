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
    #: n x 2, in the order of `dims`
    vertices: np.ndarray
    inside: bool = True

    def contains(self, column: Column) -> np.ndarray:
        mask = points_in_polygon(column(self.dims[0]), column(self.dims[1]),
                                 np.asarray(self.vertices, dtype=np.float64))
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
