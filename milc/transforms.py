"""
Axis transforms: how a measured value becomes a position on a plot, and
what the algorithms are fed.

`Biex` follows FlowJo's biexponential through the routine cytolib
translated from it. It has NOT been checked here against values out of
FlowJo itself — only its shape (monotonic, zero where expected, full range
at the top). Gates read from a workspace do not depend on it unless a
polygon is drawn on a biexponential axis.

`Logicle` is the published transform (Parks, Roederer & Moore 2006, in
the T, W, M, A of Moore & Parks 2012), checked against FlowUtils' C
translation of the authors' own code: see `test_logicle_reference`.
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


class Logicle:
    """The Logicle scale of Parks, Roederer & Moore (Cytometry A 2006), in
    the parameters of Moore & Parks (Cytometry A 2012): `T`, the value at
    the top of the scale; `M`, the decades the scale would span if it were
    a plain logarithm; `W`, the decades given to the nearly linear part
    around zero; `A`, extra negative decades. 0 to 1, like Gating-ML's.

    The inverse is the closed form `B(y) = a·e^(by) − c·e^(−dy) + f`
    above the centre of the linear part, mirrored below it. The forward
    has no closed form: a table gives a first guess and Newton's method
    finishes it, to 1e-12 of the scale or better."""

    def __init__(self, T: float = 262144.0, W: float = 0.5, M: float = 4.5,
                 A: float = 0.0):
        self.T, self.W, self.M, self.A = float(T), float(W), float(M), float(A)
        if not (self.T > 0 and self.M > 0):
            raise ValueError("a Logicle scale needs T and M above zero")
        if not 0 <= self.W <= self.M / 2:
            raise ValueError("the W of a Logicle scale is between 0 and M / 2")
        if not -self.W <= self.A <= self.M - 2 * self.W:
            raise ValueError("the A of a Logicle scale is between −W and M − 2W")
        w = self.W / (self.M + self.A)
        x2 = self.A / (self.M + self.A)
        self._x1 = x2 + w            # where the data is zero
        x0 = x2 + 2 * w
        self._b = b = (self.M + self.A) * np.log(10.0)
        self._d = d = _log_root(b, w)
        c_a = np.exp(x0 * (b + d))
        mf_a = np.exp(b * self._x1) - c_a / np.exp(d * self._x1)
        self._a = self.T / ((np.exp(b) - mf_a) - c_a / np.exp(d))
        self._c = c_a * self._a
        self._f = -mf_a * self._a
        # A first guess for the forward, to start Newton from. Against
        # u = ln(1 + x / k) the curve is close to a straight line of slope
        # 1 / b at both ends, so a small table read by index (no search)
        # is already good to about 1e-8
        a, c = self._a, self._c
        self._k = (a * b * np.exp(b * self._x1) + c * d * np.exp(-d * self._x1)) / b
        ys = np.linspace(self._x1, max(1.0, self._x1) + 0.25, 32769)
        xs = self._upper(ys)
        n = 4097
        self._per_u = (n - 1) / np.log1p(xs[-1] / self._k)
        at = self._k * np.expm1(np.arange(n) / self._per_u)
        self._table = self._newton(np.interp(at, xs, ys), at)
        self._rise = np.diff(self._table)

    def _upper(self, y):
        return self._a * np.exp(self._b * y) - self._c * np.exp(-self._d * y) + self._f

    def _newton(self, y, size):
        a, b, c, d = self._a, self._b, self._c, self._d
        target = size - self._f
        for _ in range(40):
            up, down = a * np.exp(b * y), c * np.exp(-d * y)
            step = (up - down - target) / (b * up + d * down)
            y -= step
            # Newton squares the error: after a step this small the next
            # would be under 1e-13
            if not np.any(np.abs(step) > 1e-7):
                break
        return y

    def inverse(self, y):
        y = np.asarray(y, dtype=np.float64)
        below = y < self._x1
        value = self._upper(np.where(below, 2 * self._x1 - y, y))
        return np.where(below, -value, value)

    def __call__(self, x):
        x = np.asarray(x, dtype=np.float64)
        size = np.atleast_1d(np.abs(x))
        u = np.log1p(size / self._k) * self._per_u
        # past the table, its last step carried on: the curve is a plain
        # exponential there
        i = np.minimum(np.nan_to_num(u, nan=0.0, posinf=0.0), len(self._rise) - 1).astype(np.intp)
        with np.errstate(invalid="ignore", over="ignore"):   # an infinite value: not a number
            y = self._newton(self._table[i] + (u - i) * self._rise[i], size)
        return np.where(x < 0, 2 * self._x1 - y.reshape(x.shape), y.reshape(x.shape))

    def __repr__(self):
        return f"Logicle(T={self.T:g}, W={self.W:g}, M={self.M:g}, A={self.A:g})"


_SUPER = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


@dataclass(frozen=True)
class Log:
    """A plain logarithmic axis from `min` to `max`, 0 to 1. Everything at
    or below `min` — zero and the negatives a compensated marker has —
    is drawn at the bottom edge: the reason biexponential scales exist."""
    min: float = 10.0
    max: float = 262144.0

    def __call__(self, x):
        x = np.maximum(np.asarray(x, dtype=np.float64), self.min)
        return np.log10(x / self.min) / np.log10(self.max / self.min)

    def inverse(self, y):
        return self.min * (self.max / self.min) ** np.asarray(y, dtype=np.float64)


def display_range(transform, max_value: float = 262144.0) -> tuple[float, float]:
    """The ends of an axis under `transform`, in the transform's own
    coordinates."""
    if isinstance(transform, (Biex, Log, Logicle)):
        return (0.0, 1.0)
    if isinstance(transform, Arcsinh):
        # room below zero for the spread a compensated negative has
        return (float(transform(-10.0 * transform.cofactor)), float(transform(max_value)))
    if isinstance(transform, Linear):
        return (0.0, 1.0)
    return (0.0, float(max_value))


def to_spec(transform):
    """A scale as plain values, for a file; None for no scale at all."""
    if isinstance(transform, Biex):
        return {"biex": [transform.length, transform.max_range, transform.neg,
                         transform.width, transform.pos]}
    if isinstance(transform, Logicle):
        return {"logicle": [transform.T, transform.W, transform.M, transform.A]}
    if isinstance(transform, Arcsinh):
        return {"arcsinh": transform.cofactor}
    if isinstance(transform, Log):
        return {"log": [transform.min, transform.max]}
    if isinstance(transform, Linear):
        return {"linear": [transform.min, transform.max]}
    return None


def from_spec(spec):
    if not isinstance(spec, dict):
        return None
    if "biex" in spec:
        return Biex(*spec["biex"])
    if "logicle" in spec:
        return Logicle(*[float(v) for v in spec["logicle"]])
    if "arcsinh" in spec:
        return Arcsinh(float(spec["arcsinh"]))
    if "log" in spec:
        return Log(*[float(v) for v in spec["log"]])
    if "linear" in spec:
        return Linear(*[float(v) for v in spec["linear"]])
    return None


def same(a, b) -> bool:
    return to_spec(a) == to_spec(b)


def decade_ticks(transform, max_value: float = 262144.0,
                 scale: float = 1.0) -> list[tuple[float, str]]:
    """Where 0, 10², 10³… fall on a transformed axis, with their labels:
    an axis in biexponential units reads as a cytometrist expects."""
    lo, top = (v * scale for v in display_range(transform, max_value))
    # a logarithmic axis has no zero on it
    ticks = [] if isinstance(transform, Log) else [(float(transform(0.0)) * scale, "0")]
    exp = 1
    while 10.0 ** exp <= max_value:
        at = float(transform(10.0 ** exp)) * scale
        # near zero a biexponential is almost linear and 10² lands on
        # top of 0: a label is kept only where there is room for it
        room = at - ticks[-1][0] > 0.05 * (top - lo) if ticks else at > lo
        if room and (exp >= 2 or isinstance(transform, Log)):
            ticks.append((at, "10" + str(exp).translate(_SUPER)))
        exp += 1
    return ticks
