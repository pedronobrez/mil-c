"""
One sample: its file, the matrix it is compensated with, its gates — and
the one function gates use to ask for a column by name.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

from . import compensation as comp
from .fcs import FcsFile, read_fcs
from .gating import GateNode


@dataclass
class Sample:
    fcs: FcsFile
    compensation: comp.Compensation | None = None
    gates: GateNode = field(default_factory=lambda: GateNode(""))
    transforms: dict[str, object] = field(default_factory=dict)
    #: free annotations: group, day, animal… what the plots colour by
    meta: dict[str, str] = field(default_factory=dict)
    _compensated: np.ndarray | None = field(default=None, repr=False)
    _masks: dict[str, np.ndarray] | None = field(default=None, repr=False)

    @classmethod
    def open(cls, path: str, **kw) -> "Sample":
        f = read_fcs(path)
        kw.setdefault("compensation", comp.from_keywords(f.keywords))
        return cls(fcs=f, **kw)

    @property
    def name(self) -> str:
        return os.path.splitext(self.fcs.name)[0]

    def __len__(self) -> int:
        return len(self.fcs)

    @property
    def compensated(self) -> np.ndarray:
        """Every column, the fluorescence ones compensated."""
        if self._compensated is None:
            self._compensated = (
                self.fcs.events if self.compensation is None
                else self.compensation.apply(self.fcs.events, self.fcs.channels))
        return self._compensated

    def column(self, name: str) -> np.ndarray:
        """`Comp-PE-A` is compensated, `PE-A` is as acquired. A file that
        already stores `FJComp-PE-A` is read as it is."""
        if name in self.fcs.channels:
            return self.fcs.events[:, self.fcs.index(name)]
        base, wants_comp = comp.strip_prefix(name)
        source = self.compensated if wants_comp else self.fcs.events
        return source[:, self.fcs.index(base)]

    def masks(self) -> dict[str, np.ndarray]:
        if self._masks is None:
            self._masks = self.gates.apply(self.column, len(self))
        return self._masks

    def mask(self, population: str = "") -> np.ndarray:
        return self.masks()[population]

    def count(self, population: str = "") -> int:
        return int(self.mask(population).sum())


def annotate(name: str, pattern: str) -> dict[str, str]:
    """Annotations out of a file name with named groups, so the grouping
    is the user's rule and not a guess: ``14Case4`` with
    ``(?P<day>\\d+)(?P<group>[A-Za-z]+)(?P<animal>\\d+)`` gives
    day 14, group Case, animal 4. No match, no annotations."""
    m = re.match(pattern, name)
    return {k: v for k, v in m.groupdict().items() if v is not None} if m else {}
