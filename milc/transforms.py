"""
Axis transforms: how a measured value becomes a position on a plot, and
what the algorithms are fed.

`Biex` follows FlowJo's biexponential through the routine cytolib
translated from it. It has NOT been checked here against values out of
FlowJo itself — only its shape (monotonic, zero where expected, full range
at the top). Gates read from a workspace do not depend on it unless a
polygon is drawn on a biexponential axis.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Linear:
    min: float = 0.0
    max: float = 262144.0

    def __call__(self, x):
        return (np.asarray(x, dtype=np.float64) - self.min) / (self.max - self.min)

    def inverse(self, y):
        return np.asarray(y, dtype=np.float64) * (self.max - self.min) + self.min


@dataclass(frozen=True)
class Arcsinh:
    """`asinh(x / cofactor)`: 150 is the usual cofactor for conventional
    fluorescence, 5 for mass cytometry."""
    cofactor: float = 150.0

    def __call__(self, x):
        return np.arcsinh(np.asarray(x, dtype=np.float64) / self.cofactor)

    def inverse(self, y):
        return np.sinh(np.asarray(y, dtype=np.float64)) * self.cofactor


def _log_root(b: float, w: float) -> float:
    """The d in (0, b) with 2 ln d + w d = 2 ln b - w b."""
    if w == 0:
        return b
    target = 2 * np.log(b) - w * b
    lo, hi = 1e-300, b
    for _ in range(200):
        mid = (lo + hi) / 2
        if 2 * np.log(mid) + w * mid < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


class Biex:
    """FlowJo's biexponential, as a workspace describes it: `length`
    channels, `maxRange` at the top, `pos` decades, `neg` extra negative
    decades and `width` (negative: the width basis)."""

    def __init__(self, length: int = 256, max_range: float = 262144.0,
                 neg: float = 0.0, width: float = -100.0,
                 pos: float = 4.418539922):
        self.length, self.max_range = int(length), float(max_range)
        self.neg, self.width, self.pos = float(neg), float(width), float(pos)
        self._data, self._chan = self._table()

    def _table(self):
        ln10 = np.log(10.0)
        channels = 4096  # resolution of the table, not of the display
        decades = self.pos
        width = np.log10(-self.width) if self.width < 0 else 0.0
        decades -= width / 2
        extra = max(self.neg, 0.0) + width / 2
        zero = int(extra * channels / (extra + decades))
        zero = min(zero, channels // 2)
        if zero > 0:
            decades = extra * channels / zero
        width /= 2 * decades
        pos_range = ln10 * decades
        minimum = self.max_range / np.exp(pos_range)
        neg_range = _log_root(pos_range, width)
        n = channels + 1
        j = np.arange(n, dtype=np.float64)
        positive = np.exp(j / n * pos_range)
        negative = np.exp(j / n * -neg_range)
        negative *= np.exp((pos_range + neg_range) * (width + extra / decades))
        s = positive[zero] - negative[zero]
        data = np.empty(n)
        data[zero:] = minimum * (positive[zero:] - negative[zero:] - s)
        m = 2 * zero - np.arange(zero)
        data[:zero] = -data[m]
        return data, j / channels  # data value → 0..1

    def __call__(self, x):
        return np.interp(np.asarray(x, dtype=np.float64), self._data, self._chan)

    def inverse(self, y):
        return np.interp(np.asarray(y, dtype=np.float64), self._chan, self._data)

    def __repr__(self):
        return (f"Biex(length={self.length}, max_range={self.max_range:g}, "
                f"neg={self.neg:g}, width={self.width:g}, pos={self.pos:g})")


_SUPER = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def decade_ticks(transform, max_value: float = 262144.0,
                 scale: float = 1.0) -> list[tuple[float, str]]:
    """Where 0, 10², 10³… fall on a transformed axis, with their labels:
    an axis in biexponential units reads as a cytometrist expects."""
    ticks = [(float(transform(0.0)) * scale, "0")]
    top = float(transform(max_value)) * scale
    exp = 2
    while 10.0 ** exp <= max_value:
        at = float(transform(10.0 ** exp)) * scale
        # near zero a biexponential is almost linear and 10² lands on
        # top of 0: a label is kept only where there is room for it
        if at - ticks[-1][0] > 0.05 * top:
            ticks.append((at, "10" + str(exp).translate(_SUPER)))
        exp += 1
    return ticks
