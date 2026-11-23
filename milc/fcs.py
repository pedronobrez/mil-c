"""
FCS 2.0 / 3.0 / 3.1 list-mode files, read and written with numpy alone.

`read_fcs(path)` returns the events as stored, scaled only where the file
itself asks for it (integer data with `$PnE` log amplification or `$PnG`
gain). Compensation is not applied here: see `compensation.py`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np


class FcsError(ValueError):
    pass


@dataclass
class FcsFile:
    path: str
    version: str
    #: events x parameters, float32 (float64 where the file stores doubles)
    events: np.ndarray
    #: `$PnN`, the detector names: what gates and matrices refer to
    channels: list[str]
    #: `$PnS`, the stain names; "" where the file gives none
    markers: list[str]
    #: `$PnR`
    ranges: list[float]
    #: every keyword, upper-cased key
    keywords: dict[str, str] = field(default_factory=dict)

    @property
    def name(self) -> str:
        return os.path.basename(self.path)

    def __len__(self) -> int:
        return self.events.shape[0]

    def index(self, channel: str) -> int:
        try:
            return self.channels.index(channel)
        except ValueError:
            raise KeyError(f"{self.name}: no channel {channel!r}") from None

    def label(self, channel: str) -> str:
        """`CD8 (FITC-A)`: the stain where there is one, the detector alone
        where there is not."""
        marker = self.markers[self.index(channel)]
        return f"{marker} ({channel})" if marker and marker != channel else channel


# --- reading -------------------------------------------------------------

def _parse_text(raw: str) -> dict[str, str]:
    if not raw:
        raise FcsError("empty TEXT segment")
    delim, body = raw[0], raw[1:]
    # the standard escapes a delimiter inside a value by doubling it
    tokens, cur, i, n = [], [], 0, len(body)
    while i < n:
        ch = body[i]
        if ch == delim:
            if i + 1 < n and body[i + 1] == delim:
                cur.append(delim)
                i += 2
                continue
            tokens.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
        i += 1
    if cur:
        tokens.append("".join(cur))
    kv = {tokens[j].strip().upper(): tokens[j + 1]
          for j in range(0, len(tokens) - 1, 2)}
    if "$PAR" not in kv:
        # files that write empty values (two delimiters in a row) break the
        # escape rule; read them the naive way instead
        parts = body.split(delim)
        kv = {parts[j].strip().upper(): parts[j + 1]
              for j in range(0, len(parts) - 1, 2)}
    if "$PAR" not in kv:
        raise FcsError("TEXT segment has no $PAR")
    return kv


def read_header(path: str) -> dict[str, str]:
    """The keywords alone, without touching the events: cheap enough to run
    over a whole folder."""
    with open(path, "rb") as f:
        return _read_text(f)[1]


def _read_text(f) -> tuple[str, dict[str, str], int, int]:
    head = f.read(58)
    if len(head) < 58 or head[:3] != b"FCS":
        raise FcsError("not an FCS file")
    version = head[:6].decode("ascii")

    def off(a, b):
        s = head[a:b].strip()
        return int(s) if s else 0

    t0, t1, d0, d1 = off(10, 18), off(18, 26), off(26, 34), off(34, 42)
    f.seek(t0)
    kv = _parse_text(f.read(t1 - t0 + 1).decode("latin-1"))
    # a second TEXT segment carries the rest when the first one overflowed
    s0, s1 = int(kv.get("$BEGINSTEXT", 0) or 0), int(kv.get("$ENDSTEXT", 0) or 0)
    if s0 and s1 > s0 and s0 != t0:
        f.seek(s0)
        try:
            extra = _parse_text(f.read(s1 - s0 + 1).decode("latin-1"))
            kv = {**extra, **kv}
        except FcsError:
            pass
    # offsets past 99,999,999 bytes do not fit the header and are 0 there
    if not d0 or not d1:
        d0, d1 = int(kv["$BEGINDATA"]), int(kv["$ENDDATA"])
    return version, kv, d0, d1


def read_fcs(path: str) -> FcsFile:
    with open(path, "rb") as f:
        version, kv, d0, d1 = _read_text(f)
        mode = kv.get("$MODE", "L").strip().upper()
        if mode != "L":
            raise FcsError(f"$MODE {mode}: only list mode is supported")
        npar = int(kv["$PAR"])
        dtype_code = kv["$DATATYPE"].strip().upper()
        order = kv.get("$BYTEORD", "1,2,3,4").strip()
        endian = "<" if order.startswith("1") else ">"
        bits = [int(kv[f"$P{i}B"]) for i in range(1, npar + 1)]
        f.seek(d0)
        nbytes = d1 - d0 + 1
        buf = f.read(nbytes)

    if dtype_code in ("F", "D"):
        dt = np.dtype(endian + ("f4" if dtype_code == "F" else "f8"))
        row = dt.itemsize * npar
        n = len(buf) // row
        events = np.frombuffer(buf, dtype=dt, count=n * npar).reshape(n, npar)
        events = events.astype(dt.newbyteorder("="), copy=True)
    elif dtype_code == "I":
        if any(b not in (8, 16, 32, 64) for b in bits):
            raise FcsError(f"integer widths {sorted(set(bits))} are not supported")
        rec = np.dtype([(f"p{i}", f"{endian}u{b // 8}") for i, b in enumerate(bits)])
        n = len(buf) // rec.itemsize
        packed = np.frombuffer(buf, dtype=rec, count=n)
        events = np.empty((n, npar), dtype=np.float64)
        for i in range(npar):
            col = packed[f"p{i}"]
            rng = int(float(kv[f"$P{i + 1}R"]))
            # bits above the range are not data
            mask = (1 << max(rng - 1, 1).bit_length()) - 1
            if mask < (1 << bits[i]) - 1:
                col = col & mask
            events[:, i] = col
        for i in range(npar):
            rng = float(kv[f"$P{i + 1}R"])
            dec, off_ = _amplification(kv.get(f"$P{i + 1}E", "0,0"))
            if dec > 0:
                events[:, i] = (off_ or 1.0) * 10 ** (dec * events[:, i] / rng)
            else:
                gain = float(kv.get(f"$P{i + 1}G", 1) or 1)
                if gain not in (0.0, 1.0):
                    events[:, i] /= gain
        events = events.astype(np.float32)
    else:
        raise FcsError(f"$DATATYPE {dtype_code} is not supported")

    declared = kv.get("$TOT")
    if declared is not None and int(declared) < events.shape[0]:
        events = events[: int(declared)]

    return FcsFile(
        path=path, version=version, events=events,
        channels=[kv[f"$P{i}N"].strip() for i in range(1, npar + 1)],
        markers=[kv.get(f"$P{i}S", "").strip() for i in range(1, npar + 1)],
        ranges=[float(kv[f"$P{i}R"]) for i in range(1, npar + 1)],
        keywords=kv,
    )


def _amplification(text: str) -> tuple[float, float]:
    try:
        a, b = (float(x) for x in text.split(","))
        return a, b
    except ValueError:
        return 0.0, 0.0


# --- writing -------------------------------------------------------------

def write_fcs(path: str, events: np.ndarray, channels: list[str],
              markers: list[str] | None = None,
              ranges: list[float] | None = None,
              keywords: dict[str, str] | None = None) -> None:
    """FCS 3.1, 32-bit floats, little-endian: what FlowJo and every other
    reader opens. Used for the concatenated file."""
    events = np.ascontiguousarray(events, dtype="<f4")
    n, npar = events.shape
    markers = markers or [""] * npar
    if ranges is None:
        ranges = [float(max(1.0, np.ceil(events[:, i].max()) if n else 1.0))
                  for i in range(npar)]
    kv: dict[str, str] = {
        "$BEGINANALYSIS": "0", "$ENDANALYSIS": "0",
        "$BEGINSTEXT": "0", "$ENDSTEXT": "0",
        "$BYTEORD": "1,2,3,4", "$DATATYPE": "F", "$MODE": "L",
        "$NEXTDATA": "0", "$PAR": str(npar), "$TOT": str(n),
    }
    for i in range(npar):
        p = f"$P{i + 1}"
        kv[p + "B"] = "32"
        kv[p + "E"] = "0,0"
        kv[p + "N"] = channels[i]
        kv[p + "R"] = repr(float(ranges[i])).rstrip("0").rstrip(".")
        if markers[i]:
            kv[p + "S"] = markers[i]
    for k, v in (keywords or {}).items():
        if v != "":
            kv.setdefault(k, str(v))

    delim = "|"

    def text(begin: int, end: int) -> bytes:
        items = {"$BEGINDATA": str(begin), "$ENDDATA": str(end), **kv}
        return (delim + delim.join(
            f"{k}{delim}{str(v).replace(delim, delim * 2)}"
            for k, v in items.items()) + delim).encode("latin-1", "replace")

    # the offsets are part of the text they sit after: settle by iterating
    t0 = 58
    begin = end = 0
    for _ in range(5):
        body = text(begin, end)
        nb = t0 + len(body)
        ne = nb + events.nbytes - 1
        if (nb, ne) == (begin, end):
            break
        begin, end = nb, ne
    body = text(begin, end)
    t1 = t0 + len(body) - 1

    def h(v: int) -> bytes:
        return (b"%8d" % v) if v <= 99_999_999 else b"       0"

    header = b"FCS3.1    " + h(t0) + h(t1) + h(begin) + h(end) + h(0) + h(0)
    with open(path, "wb") as f:
        f.write(header)
        f.write(body)
        f.write(events.tobytes())
