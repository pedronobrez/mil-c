"""
Spillover matrices: from the file's `$SPILL`, or from a FlowJo workspace
(which is the one to trust when the user edited it there).

Row i of a spillover matrix is how much fluorochrome i shows in every
detector, so `observed = true @ S` and `true = observed @ inv(S)`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: how compensated channels are named in gates. FlowJo writes `Comp-` in a
#: workspace and `FJComp-` in the files it exports
PREFIXES = ("Comp-", "FJComp-")


@dataclass
class Compensation:
    channels: list[str]
    #: len(channels) square, row = fluorochrome, column = detector
    matrix: np.ndarray
    name: str = ""

    def apply(self, events: np.ndarray, all_channels: list[str]) -> np.ndarray:
        """A copy of `events` with the matrix's channels compensated and
        every other column untouched."""
        idx = [all_channels.index(c) for c in self.channels]
        out = np.array(events, dtype=np.float32, copy=True)
        inv = np.linalg.inv(self.matrix.astype(np.float64))
        out[:, idx] = (events[:, idx].astype(np.float64) @ inv).astype(np.float32)
        return out


def from_keywords(keywords: dict[str, str]) -> Compensation | None:
    for key in ("$SPILLOVER", "SPILL", "$SPILL", "SPILLOVER"):
        text = keywords.get(key)
        if text:
            break
    else:
        return None
    parts = [p.strip() for p in text.split(",")]
    n = int(parts[0])
    names = parts[1:1 + n]
    values = np.array([float(v) for v in parts[1 + n:1 + n + n * n]])
    if values.size != n * n:
        return None
    return Compensation(names, values.reshape(n, n), name="acquisition")


def strip_prefix(channel: str) -> tuple[str, bool]:
    """`Comp-PE-A` → (`PE-A`, True); `FSC-A` → (`FSC-A`, False)."""
    for p in PREFIXES:
        if channel.startswith(p):
            return channel[len(p):], True
    return channel, False
