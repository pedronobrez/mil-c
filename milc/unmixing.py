"""
Spectral unmixing: for cytometers with many detectors and fewer
fluorochromes (Cytek Aurora and Northern Lights, Sony ID7000, BD
FACSymphony A5 SE).

Every fluorochrome has a spectrum — how much of it each detector sees,
highest detector = 1. An event's detector values `y` are taken to be a
sum of the spectra, `y ≈ a @ M`, and unmixing finds the abundances `a`:

* ordinary least squares, `a = y @ pinv(M)`: one matrix for every event;
* weighted least squares, as Cytek's WLS: each detector of each event
  weighted by 1 / its variance, taken as Poisson-like — the value
  itself, never below a floor — so that a bright detector, which is
  noisy, counts for less than a dim one. One small system per event.

The spectra come from single-stain reference controls (the positive
events' median minus the negative's — the same control's dim events or
an unstained control —, normalised), optionally with autofluorescence as
one more spectrum, from an unstained control of the cells; or from a CSV.

The unmixed channels are named after the fluorochromes and replace the
detectors in the sample (`apply_fcs`): scatter, time and anything else
the spectra do not cover stay as they were. To the rest of the program
they are compensated fluorescence, as though by a matrix: see
`Unmixing.identity` and `Session.set_unmixing`.

How well the panel can be unmixed is said by the similarity matrix
(cosine similarity of every two spectra, as Cytek's) and the complexity
index (the condition number of the matrix of spectra, the largest
singular value over the smallest — AutoSpectral's
`calculate.condition.number`; Cytek's own formula is not published).
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field

import numpy as np

from .compensation import Compensation
from .fcs import FcsFile

#: AutoSpectral's default floor for the weights (raw Aurora units): below
#: this, a detector's variance is taken as the floor, not as its value
NOISE_FLOOR = 125.0
OLS, WLS = "ols", "wls"
AF = "AF"

_NOT_DETECTOR = re.compile(r"^(FSC|SSC|TIME|EVENT|WIDTH|HEIGHT|SAMPLEID|ORIGINAL|"
                           r"CLUSTER|UMAP|TSNE|FJCOMP|COMP-)", re.IGNORECASE)
_SUFFIX = re.compile(r"-(A|H|W)$")


def detector_channels(channels: list[str]) -> list[str]:
    """The fluorescence detectors of a spectral file, in its order: every
    channel but scatter, time and the like; of those, the areas (`-A`)
    alone when the file has areas."""
    found = [c for c in channels if not _NOT_DETECTOR.match(c.strip())]
    areas = [c for c in found if c.endswith("-A")]
    return areas or found


def _base(channel: str) -> str:
    return _SUFFIX.sub("", channel)


# --- the spectra of reference controls ----------------------------------------


@dataclass
class Reference:
    """One spectrum and where it came from: what the dialog shows."""
    fluor: str
    spectrum: np.ndarray
    marker: str = ""
    source: str = ""
    #: the detector the spectrum peaks on
    peak: str = ""
    positive_events: int = 0
    negative_events: int = 0
    #: median of the positive events on the peak detector, and of the
    #: negative ones: how far apart they are
    positive_median: float = 0.0
    negative_median: float = 0.0


def _usable_rows(events: np.ndarray, ranges: list[float] | None) -> np.ndarray:
    """Events not at the top of any detector's range: a saturated detector
    bends the spectrum."""
    if ranges is None:
        return np.ones(events.shape[0], dtype=bool)
    top = np.asarray(ranges, dtype=float) * 0.99
    return ~(events >= top).any(axis=1)


def spectrum_of(events: np.ndarray, negative: np.ndarray | None = None,
                ranges: list[float] | None = None, positive_fraction: float = 0.05,
                negative_fraction: float = 0.3) -> dict:
    """The spectrum of one single-stain control.

    `events` is events x detectors. The peak detector is the one where the
    brightest events stand furthest above the negative; the positive
    events are the brightest `positive_fraction` on it (saturated ones
    left out), the negative ones the dimmest `negative_fraction` of the
    same control — or `negative`, events of an unstained control, when
    given. The spectrum is the positive median minus the negative median,
    divided by its largest value."""
    events = np.asarray(events, dtype=np.float64)
    ok = _usable_rows(events, ranges)
    x = events[ok]
    if x.shape[0] < 20:
        raise ValueError("fewer than 20 events below saturation")
    if negative is not None:
        neg = np.asarray(negative, dtype=np.float64)
        neg = neg[_usable_rows(neg, ranges)]
        neg_median = np.median(neg, axis=0)
    else:
        neg_median = np.median(x, axis=0)    # only to choose the peak
    bright = np.quantile(x, 0.995, axis=0)
    peak = int(np.argmax(bright - neg_median))
    on_peak = x[:, peak]
    pos_rows = on_peak >= np.quantile(on_peak, 1.0 - positive_fraction)
    if negative is None:
        neg = x[on_peak <= np.quantile(on_peak, negative_fraction)]
        neg_median = np.median(neg, axis=0)
    pos_median = np.median(x[pos_rows], axis=0)
    raw = pos_median - neg_median
    top = raw.max()
    if not np.isfinite(top) or top <= 0:
        raise ValueError("the positive events are not brighter than the negative ones")
    return {"spectrum": raw / top, "peak": peak, "positive_events": int(pos_rows.sum()),
            "negative_events": int(neg.shape[0]), "positive_median": float(pos_median[peak]),
            "negative_median": float(neg_median[peak])}


def autofluorescence_of(events: np.ndarray, ranges: list[float] | None = None) -> np.ndarray:
    """The autofluorescence spectrum of unstained cells: their median,
    divided by its largest value."""
    events = np.asarray(events, dtype=np.float64)
    m = np.median(events[_usable_rows(events, ranges)], axis=0)
    if m.max() <= 0:
        raise ValueError("the unstained control has no signal above zero")
    return m / m.max()


def noise_of(events: np.ndarray) -> np.ndarray:
    """Per detector, a floor for the WLS weights: the variance of
    unstained events, robustly (MAD²), and never below `NOISE_FLOOR`."""
    events = np.asarray(events, dtype=np.float64)
    mad = 1.4826 * np.median(np.abs(events - np.median(events, axis=0)), axis=0)
    return np.maximum(mad, NOISE_FLOOR)


def read_references(controls: list[dict], detectors: list[str] | None = None,
                    unstained: str | None = None, autofluorescence: bool = False,
                    positive_fraction: float = 0.05, negative_fraction: float = 0.3,
                    reader=None) -> tuple[list[str], list[Reference], np.ndarray | None]:
    """Spectra from files. `controls` is [{"path", "fluor", "marker"?,
    "negative"?}], `negative` being "internal" (the control's own dim
    events, the default) or "unstained" (the `unstained` file's events).
    With `autofluorescence`, the unstained control's own spectrum is one
    more, named `AF`. Returns (detectors, references, noise floor or None)."""
    from .fcs import read_fcs
    reader = reader or read_fcs
    cache: dict[str, FcsFile] = {}

    def load(path):
        if path not in cache:
            cache[path] = reader(path)
        return cache[path]

    if not controls:
        raise ValueError("no reference controls")
    first = load(controls[0]["path"])
    detectors = list(detectors or detector_channels(first.channels))

    def block(f: FcsFile):
        missing = [d for d in detectors if d not in f.channels]
        if missing:
            raise ValueError(f"{f.name} has no {', '.join(missing[:5])}")
        idx = [f.index(d) for d in detectors]
        return f.events[:, idx], [f.ranges[i] for i in idx]

    neg_events, noise = None, None
    if unstained:
        neg_events, _ = block(load(unstained))
        noise = noise_of(neg_events)
    refs = []
    for c in controls:
        f = load(c["path"])
        events, ranges = block(f)
        use_unstained = c.get("negative", "internal") == "unstained"
        if use_unstained and neg_events is None:
            raise ValueError(f"{c['fluor']}: an unstained control is needed for its negative")
        found = spectrum_of(events, neg_events if use_unstained else None, ranges,
                            positive_fraction, negative_fraction)
        refs.append(Reference(c["fluor"], found["spectrum"], c.get("marker", ""), f.name,
                              detectors[found["peak"]], found["positive_events"],
                              found["negative_events"], found["positive_median"],
                              found["negative_median"]))
    if autofluorescence:
        if unstained is None:
            raise ValueError("autofluorescence needs an unstained control")
        f = load(unstained)
        events, ranges = block(f)
        af = autofluorescence_of(events, ranges)
        refs.append(Reference(AF, af, "Autofluorescence", f.name,
                              detectors[int(np.argmax(af))], int(events.shape[0])))
    names = [r.fluor for r in refs]
    if len(set(names)) != len(names):
        raise ValueError("two references have the same name")
    return detectors, refs, noise


# --- how well the panel unmixes -----------------------------------------------


def similarity_matrix(spectra: np.ndarray) -> np.ndarray:
    """Cosine similarity of every two spectra (rows): 1 is the same shape,
    0 nothing in common. Cytek warns above 0.98."""
    s = np.asarray(spectra, dtype=np.float64)
    norm = np.sqrt((s ** 2).sum(axis=1))
    return (s @ s.T) / np.outer(norm, norm)


def complexity_index(spectra: np.ndarray) -> float:
    """The condition number of the matrix of spectra, each with its peak at
    1: the largest singular value over the smallest. The more alike some
    combination of the spectra, the larger it is, and the more each
    fluorochrome's noise spreads into the others."""
    s = np.asarray(spectra, dtype=np.float64)
    s = s / s.max(axis=1, keepdims=True)
    sv = np.linalg.svd(s, compute_uv=False)
    return float(sv.max() / sv.min()) if sv.min() > 0 else float("inf")


# --- unmixing -------------------------------------------------------------------


@dataclass
class Unmixing:
    detectors: list[str]
    #: the unmixed channels' names, one per spectrum
    fluors: list[str]
    #: len(fluors) x len(detectors), each row's largest value 1
    spectra: np.ndarray
    #: the stain on each fluorochrome; "" where none was given
    markers: list[str] = field(default_factory=list)
    method: str = OLS
    #: WLS: per detector, the variance below which a value is not trusted
    noise: np.ndarray | None = None
    #: the file each spectrum came from, to say so
    sources: list[str] = field(default_factory=list)

    def __post_init__(self):
        self.spectra = np.asarray(self.spectra, dtype=np.float64)
        k, d = self.spectra.shape
        if k != len(self.fluors) or d != len(self.detectors):
            raise ValueError(f"{k} x {d} spectra for {len(self.fluors)} fluorochromes "
                             f"and {len(self.detectors)} detectors")
        if k > d:
            raise ValueError(f"{k} fluorochromes cannot be unmixed from {d} detectors")
        if np.linalg.matrix_rank(self.spectra) < k:
            raise ValueError("some spectra are a combination of the others")
        if len(set(self.fluors)) != k:
            raise ValueError("two fluorochromes have the same name")
        self.markers = list(self.markers) + [""] * (k - len(self.markers))
        if self.method not in (OLS, WLS):
            raise ValueError(f"unknown method {self.method!r}")

    @classmethod
    def from_references(cls, detectors: list[str], refs: list[Reference],
                        method: str = OLS, noise=None) -> "Unmixing":
        return cls(list(detectors), [r.fluor for r in refs],
                   np.array([r.spectrum for r in refs]), [r.marker for r in refs],
                   method, None if noise is None else np.asarray(noise, dtype=float),
                   [r.source for r in refs])

    # --- the numbers -------------------------------------------------------

    def ols_matrix(self) -> np.ndarray:
        """detectors x fluorochromes: `a = y @ U`."""
        return np.linalg.pinv(self.spectra)

    def unmix(self, y: np.ndarray, method: str | None = None,
              chunk: int = 4096) -> np.ndarray:
        """events x detectors (in `detectors` order) → events x fluorochromes."""
        y = np.asarray(y, dtype=np.float64)
        if (method or self.method) == OLS:
            return y @ self.ols_matrix()
        floor = NOISE_FLOOR if self.noise is None else np.asarray(self.noise, dtype=float)
        M = self.spectra
        k = M.shape[0]
        # every event's M W Mᵀ at once, as one product: entry (i, j) is
        # Σ_d w_d M_id M_jd, so the weights times the products of rows
        iu = np.triu_indices(k)
        pairs = (M[iu[0]] * M[iu[1]]).T                        # d x k(k+1)/2
        out = np.empty((y.shape[0], k))
        for start in range(0, y.shape[0], chunk):
            yy = y[start:start + chunk]
            w = 1.0 / np.maximum(np.abs(yy), floor)          # n x d
            gram = np.empty((yy.shape[0], k, k))
            upper = w @ pairs
            gram[:, iu[0], iu[1]] = upper
            gram[:, iu[1], iu[0]] = upper
            rhs = (w * yy) @ M.T                               # n x k
            out[start:start + chunk] = np.linalg.solve(gram, rhs[:, :, None])[:, :, 0]
        return out

    def residual(self, y: np.ndarray, a: np.ndarray | None = None) -> np.ndarray:
        """Per event, how much of it the spectra do not explain:
        |y − a M| / |y|."""
        y = np.asarray(y, dtype=np.float64)
        a = self.unmix(y) if a is None else a
        left = np.linalg.norm(y - a @ self.spectra, axis=1)
        return left / np.maximum(np.linalg.norm(y, axis=1), 1e-9)

    def similarity(self) -> np.ndarray:
        return similarity_matrix(self.spectra)

    def complexity(self) -> float:
        return complexity_index(self.spectra)

    # --- the program ---------------------------------------------------------

    def covers(self, channels: list[str]) -> bool:
        return all(d in channels for d in self.detectors)

    def kept(self, channels: list[str]) -> list[str]:
        """The channels that stay as they are: neither a detector of the
        spectra nor its height or width."""
        gone = {_base(d) for d in self.detectors}
        return [c for c in channels if _base(c) not in gone]

    def channels_of(self, channels: list[str]) -> list[str]:
        return self.kept(channels) + list(self.fluors)

    def apply_fcs(self, f: FcsFile) -> FcsFile:
        """The file as the program sees it once unmixed: the kept channels
        as acquired, then one channel per fluorochrome."""
        missing = [d for d in self.detectors if d not in f.channels]
        if missing:
            raise KeyError(f"{f.name}: no {', '.join(missing[:5])}: "
                           "the spectra are not for this file")
        kept = self.kept(f.channels)
        ki = [f.index(c) for c in kept]
        di = [f.index(d) for d in self.detectors]
        a = self.unmix(f.events[:, di]).astype(np.float32)
        events = np.column_stack([f.events[:, ki].astype(np.float32), a])
        top = max(f.ranges[i] for i in di)
        keywords = {k: v for k, v in f.keywords.items()
                    if k not in ("$SPILLOVER", "SPILL", "$SPILL", "SPILLOVER")}
        keywords["MILC UNMIXED"] = self.method.upper()
        return FcsFile(path=f.path, version=f.version, events=events,
                       channels=kept + list(self.fluors),
                       markers=[f.markers[i] for i in ki]
                       + [m or "" for m in self.markers],
                       ranges=[f.ranges[i] for i in ki] + [top] * len(self.fluors),
                       keywords=keywords)

    def identity(self) -> Compensation:
        """What the unmixed channels are compensated with: nothing more —
        but saying so makes them fluorescence to everything else (the
        biexponential axis, `Comp-` in gates, normalisation, statistics)."""
        return Compensation(list(self.fluors), np.eye(len(self.fluors)), "unmixed")

    # --- saved with the project ----------------------------------------------

    def to_dict(self) -> dict:
        return {"detectors": self.detectors, "fluors": self.fluors,
                "spectra": self.spectra.tolist(), "markers": self.markers,
                "method": self.method,
                "noise": None if self.noise is None else np.asarray(self.noise).tolist(),
                "sources": self.sources}

    @classmethod
    def from_dict(cls, data) -> "Unmixing | None":
        if not data:
            return None
        try:
            return cls(list(data["detectors"]), list(data["fluors"]),
                       np.array(data["spectra"], dtype=float), list(data.get("markers") or []),
                       str(data.get("method") or OLS),
                       None if data.get("noise") is None else np.array(data["noise"], dtype=float),
                       list(data.get("sources") or []))
        except (KeyError, TypeError, ValueError):
            return None


# --- spectra in a CSV -------------------------------------------------------------


def read_csv(path: str, detectors: list[str] | None = None
             ) -> tuple[list[str], list[str], np.ndarray]:
    """(detectors, fluorochromes, spectra) from a table with one
    fluorochrome per row and one detector per column (a header row, the
    names in the first column) — or the other way round, which is told
    from where the detector names are. Each spectrum is normalised to a
    peak of 1. With `detectors`, the table is put in their order."""
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = [r for r in csv.reader(f) if r and any(x.strip() for x in r)]
    if len(rows) < 2:
        raise ValueError("the table has no rows of numbers")
    head = [h.strip() for h in rows[0][1:]]
    names = [r[0].strip() for r in rows[1:]]
    values = np.array([[float(x) if x.strip() else 0.0 for x in r[1:1 + len(head)]]
                       for r in rows[1:]], dtype=float)
    looks_like_detectors = (lambda xs: sum(bool(re.match(r"^[A-Za-z]+\d+(-[AHW])?$", x))
                                           for x in xs) > len(xs) / 2)
    by_column = (detectors is not None and sum(d in head for d in detectors)
                 >= sum(d in names for d in detectors)) or \
        (detectors is None and looks_like_detectors(head) >= looks_like_detectors(names))
    if by_column:
        dets, fluors, spectra = head, names, values
    else:
        dets, fluors, spectra = names, head, values.T
    if detectors is not None:
        missing = [d for d in detectors if d not in dets]
        if missing:
            raise ValueError(f"the table has no {', '.join(missing[:5])}")
        spectra = spectra[:, [dets.index(d) for d in detectors]]
        dets = list(detectors)
    top = spectra.max(axis=1, keepdims=True)
    if (top <= 0).any():
        raise ValueError("a spectrum has no value above zero")
    return list(dets), list(fluors), spectra / top


def write_csv(path: str, unmixing: Unmixing) -> None:
    """One fluorochrome per row, one detector per column: what `read_csv`
    and spreadsheets read."""
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Fluorochrome", *unmixing.detectors])
        for name, row in zip(unmixing.fluors, unmixing.spectra):
            w.writerow([name, *(f"{v:.6g}" for v in row)])
