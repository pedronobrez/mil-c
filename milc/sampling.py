"""
How many events each sample contributes to the concatenated file, decided
before anything is drawn.

A treated sample usually carries several times the cells of its control, and
an embedding built on everything is an embedding of the treated group.
The plan says what each sample gives, which sample set the limit, and what
was left out — so the cost of balancing is read before it is paid.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass, field

import numpy as np

#: the same number from every sample: the smallest one sets it
EQUAL_PER_SAMPLE = "equal_per_sample"
#: the same total from every group, spread as evenly as its samples allow
EQUAL_PER_GROUP = "equal_per_group"
#: at most n from each; samples below n give all they have
CEILING = "ceiling"
MODES = (EQUAL_PER_SAMPLE, EQUAL_PER_GROUP, CEILING)


@dataclass
class Plan:
    mode: str
    #: sample → events to draw (0 for an excluded sample)
    take: dict[str, int]
    available: dict[str, int]
    groups: dict[str, str]
    #: sample → why it gives nothing
    excluded: dict[str, str] = field(default_factory=dict)
    #: the sample (or group) whose size set everyone's number
    limited_by: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def total(self) -> int:
        return sum(self.take.values())

    def group_totals(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for s, k in self.take.items():
            g = self.groups.get(s, "")
            out[g] = out.get(g, 0) + k
        return out

    def used_fraction(self) -> float:
        """Of the events the included samples have, the share kept."""
        have = sum(a for s, a in self.available.items() if s not in self.excluded)
        return self.total / have if have else 0.0


def plan(available: dict[str, int], groups: dict[str, str] | None = None,
         mode: str = EQUAL_PER_SAMPLE, n: int | None = None,
         min_events: int = 1) -> Plan:
    """`n` is per sample for the first and third modes and per group for
    the second; left out, the data sets it. A sample with fewer than
    `min_events` is excluded rather than allowed to drag everyone down."""
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r}")
    groups = dict(groups or {})
    min_events = max(int(min_events), 1)
    excluded = {s: f"{a} events, below the minimum of {min_events}"
                for s, a in available.items() if a < min_events}
    ok = {s: a for s, a in available.items() if s not in excluded}
    take = {s: 0 for s in available}
    p = Plan(mode, take, dict(available), groups, excluded)
    if not ok:
        p.notes.append("no sample has enough events")
        return p

    if mode == EQUAL_PER_SAMPLE:
        smallest = min(ok, key=lambda s: (ok[s], s))
        k = ok[smallest]
        if n is not None and n < k:
            k = n
        else:
            p.limited_by = smallest
            if n is not None and n > k:
                p.notes.append(
                    f"asked for {n} per sample; {smallest} has only {k}")
        for s in ok:
            take[s] = k

    elif mode == CEILING:
        if n is None:
            raise ValueError("the ceiling mode needs n")
        for s, a in ok.items():
            take[s] = min(a, n)
        short = [s for s, a in ok.items() if a < n]
        if short:
            p.notes.append(
                f"{len(short)} sample(s) have fewer than {n} and give all "
                "they have: the samples are not balanced")

    else:  # EQUAL_PER_GROUP
        members: dict[str, list[str]] = {}
        for s in ok:
            members.setdefault(groups.get(s, ""), []).append(s)
        sums = {g: sum(ok[s] for s in ss) for g, ss in members.items()}
        smallest = min(sums, key=lambda g: (sums[g], g))
        target = sums[smallest]
        if n is not None and n < target:
            target = n
        else:
            p.limited_by = smallest
        for g, ss in members.items():
            for s, k in _water_fill({s: ok[s] for s in ss}, target).items():
                take[s] = k
        if any(len({take[s] for s in ss}) > 1 for ss in members.values()):
            p.notes.append(
                "inside a group the samples do not give the same number: "
                "the larger ones make up for the smaller")
    return p


def _water_fill(have: dict[str, int], target: int) -> dict[str, int]:
    """`target` events out of these samples, as evenly as they allow: the
    small ones give everything, the rest share what is left."""
    out = {s: 0 for s in have}
    left = min(target, sum(have.values()))
    pending = sorted(have, key=lambda s: (have[s], s))
    while pending:
        share = left // len(pending)
        s = pending[0]
        if have[s] <= share:
            out[s] = have[s]
            left -= have[s]
            pending.pop(0)
            continue
        extra = left - share * len(pending)
        # the odd events go to the largest samples, always the same ones
        for i, s in enumerate(reversed(pending)):
            out[s] = share + (1 if i < extra else 0)
        break
    return out


def draw(mask: np.ndarray, k: int, seed: int, key: str) -> np.ndarray:
    """`k` event indices out of those `mask` selects, in acquisition
    order. The draw depends on the seed and the sample's own name, not on
    which other samples are in the batch or in what order."""
    pool = np.flatnonzero(mask)
    if k >= pool.size:
        return pool
    rng = np.random.default_rng([int(seed), zlib.crc32(key.encode("utf-8"))])
    return np.sort(rng.choice(pool, size=k, replace=False))
