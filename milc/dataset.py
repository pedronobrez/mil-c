"""
The concatenated events every analysis runs on, with where each one came
from: which sample, which event of that sample.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from .compensation import strip_prefix
from .fcs import read_fcs, write_fcs
from .sample import Sample
from .sampling import Plan, draw


@dataclass
class Dataset:
    #: events x channels. Fluorescence channels are compensated, nothing
    #: is transformed: transforms are chosen where the data is used
    X: np.ndarray
    channels: list[str]
    markers: list[str]
    #: which channels the compensation touched
    compensated: list[bool]
    #: per event: index into `samples`
    sample_index: np.ndarray
    #: per event: its position in the original file, from 0
    event_index: np.ndarray
    samples: list[str]
    #: per sample: group, day, animal…
    meta: list[dict[str, str]]
    #: how this was built: population, sampling plan, seed
    provenance: dict = field(default_factory=dict)
    #: per-event results by name: `UMAP1`, `UMAP2`, `FlowSOM`…
    derived: dict[str, np.ndarray] = field(default_factory=dict)
    #: how each result was computed, by the same names
    runs: dict[str, dict] = field(default_factory=dict)

    def __len__(self) -> int:
        return self.X.shape[0]

    def column(self, channel: str) -> np.ndarray:
        return self.X[:, self.channels.index(channel)]

    def obs(self, key: str) -> np.ndarray:
        """A per-event annotation: `sample`, or any key of the metadata."""
        values = self.samples if key == "sample" else [m.get(key, "") for m in self.meta]
        return np.asarray(values, dtype=object)[self.sample_index]

    def matrix(self, channels: list[str], transform=None,
               scale: float = 1.0) -> np.ndarray:
        """What an algorithm is fed: these channels, transformed. `transform`
        is one callable for all of them or a dict by channel; FlowJo hands
        its plugins `Biex()` on a scale of 1023 (measured: see
        tests/test_real.py)."""
        cols = []
        for c in channels:
            t = transform.get(c) if isinstance(transform, dict) else transform
            x = self.column(c).astype(np.float64)
            cols.append(t(x) * scale if t is not None else x * scale)
        return np.column_stack(cols)

    def resolve(self, name: str) -> str:
        """A channel from what a person types: the detector (`FITC-A`,
        with or without `Comp-`), or the stain or its first word (`CD8`),
        whatever the case. Ambiguity is an error, not a guess."""
        want = strip_prefix(name.strip())[0].lower()
        for c in self.channels:
            if c.lower() == want:
                return c
        hits = [c for c, m in zip(self.channels, self.markers)
                if m and (m.lower() == want or m.lower().split()[0] == want)]
        if len(hits) != 1:
            raise KeyError(f"{name!r}: " + ("matches " + ", ".join(hits) if hits
                                            else "no such channel or stain"))
        return hits[0]

    def label(self, channel: str) -> str:
        marker = self.markers[self.channels.index(channel)]
        return marker or channel

    def counts(self) -> dict[str, int]:
        n = np.bincount(self.sample_index, minlength=len(self.samples))
        return dict(zip(self.samples, (int(v) for v in n)))

    # --- files ----------------------------------------------------------

    def save(self, path: str) -> None:
        np.savez_compressed(
            path, X=self.X, sample_index=self.sample_index,
            event_index=self.event_index,
            **{"d:" + k: v for k, v in self.derived.items()},
            header=json.dumps({
                "runs": self.runs,
                "channels": self.channels, "markers": self.markers,
                "compensated": self.compensated, "samples": self.samples,
                "meta": self.meta, "provenance": self.provenance}))

    @classmethod
    def load(cls, path: str) -> "Dataset":
        with np.load(path, allow_pickle=False) as z:
            h = json.loads(str(z["header"]))
            derived = {k[2:]: z[k] for k in z.files if k.startswith("d:")}
            return cls(X=z["X"], sample_index=z["sample_index"],
                       event_index=z["event_index"], derived=derived, **h)

    @classmethod
    def from_fcs(cls, path: str, sample_column: str | None = "SampleID") -> "Dataset":
        """A file that is already concatenated — FlowJo's, or our own
        export. `FJComp-`/`Comp-` channels are taken as compensated and
        lose the prefix; `sample_column`, where present, says which
        sample each event came from (any distinct values, in order)."""
        f = read_fcs(path)
        keep = [i for i, c in enumerate(f.channels) if c != sample_column]
        stripped = [strip_prefix(f.channels[i]) for i in keep]
        if sample_column in f.channels:
            ids, sidx = np.unique(f.events[:, f.index(sample_column)],
                                  return_inverse=True)
            named = f.keywords.get("MILC SAMPLES", "").split(";")
            samples = named if len(named) == ids.size else [f"{v:g}" for v in ids]
        else:
            sidx, samples = np.zeros(len(f), dtype=np.int64), [f.name]
        return cls(
            X=f.events[:, keep].astype(np.float32),
            channels=[c for c, _ in stripped],
            markers=[f.markers[i] for i in keep],
            compensated=[flag for _, flag in stripped],
            sample_index=sidx.astype(np.int32),
            event_index=np.arange(len(f), dtype=np.int64),
            samples=samples, meta=[{} for _ in samples],
            provenance={"source": f.name})

    def to_fcs(self, path: str) -> None:
        """For FlowJo or anything else: compensated channels are named
        `Comp-…` and `SampleID` (from 1, in the order of `samples`) is the
        last column."""
        names = [("Comp-" + c) if comp else c
                 for c, comp in zip(self.channels, self.compensated)]
        data = np.column_stack([self.X, self.sample_index.astype(np.float32) + 1])
        ranges = [262144.0] * len(names) + [float(len(self.samples) + 1)]
        keywords = {"MILC SAMPLES": ";".join(self.samples)}
        cleaned = self.provenance.get("quality")
        if cleaned:
            # events removed by quality control / events of the file, per sample
            lost = cleaned["samples"]
            keywords["MILC QC"] = ";".join(
                f"{lost[name]['removed']}/{lost[name]['events']}" if name in lost else ""
                for name in self.samples)
        write_fcs(path, data, names + ["SampleID"], self.markers + [""], ranges, keywords)


def concatenate(samples: list[Sample], plan: Plan, population: str = "",
                seed: int = 0) -> Dataset:
    """Draw what the plan says from `population` of every sample and stack
    it. Samples must share their channels, in order."""
    used = [s for s in samples if plan.take.get(s.name, 0) > 0]
    if not used:
        raise ValueError("the plan takes no events from any sample")
    first = used[0]
    for s in used[1:]:
        if s.fcs.channels != first.fcs.channels:
            raise ValueError(f"{s.name} and {first.name} do not have the same channels")
    blocks, sidx, eidx = [], [], []
    for i, s in enumerate(used):
        rows = draw(s.mask(population), plan.take[s.name], seed, s.name)
        blocks.append(s.compensated[rows])
        sidx.append(np.full(rows.size, i, dtype=np.int32))
        eidx.append(rows.astype(np.int64))
    touched = set(first.compensation.channels) if first.compensation else set()
    return Dataset(
        X=np.concatenate(blocks).astype(np.float32),
        channels=list(first.fcs.channels), markers=list(first.fcs.markers),
        compensated=[c in touched for c in first.fcs.channels],
        sample_index=np.concatenate(sidx), event_index=np.concatenate(eidx),
        samples=[s.name for s in used], meta=[dict(s.meta) for s in used],
        provenance={
            "population": population, "seed": int(seed), "mode": plan.mode,
            "take": {s.name: plan.take[s.name] for s in used},
            "available": plan.available, "excluded": plan.excluded,
            "limited_by": plan.limited_by,
        })
