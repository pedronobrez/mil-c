"""What the four comparisons share: where the work goes, how R is called,
how the real experiment is read. Nothing here is collected by pytest."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))     # the repository: `milc`

#: everything written goes here, never into the repository
WORK = os.environ.get("MILC_REF_WORK") or os.path.join(tempfile.gettempdir(), "milc_reference")
#: a folder with the experiment's .fcs files and the FlowJo workspace; the
#: comparisons on real data are left out when it is not given
DATA = os.environ.get("MILC_REF_DATA", "")
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")

PATTERN = r"(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)"
LIVE = "Cells/Single Cells/Single Cells/Live"
CD45 = LIVE + "/CD45"
#: the stains that tell cell types apart inside CD45+ (not the viability
#: dye, not CD45 itself)
MARKERS = ["FITC-A", "PE-A", "PerCP-Cy5.5-A", "APC-A", "BV421-A", "BV510-A"]


def work(name: str) -> str:
    path = os.path.join(WORK, name)
    os.makedirs(path, exist_ok=True)
    return path


def rscript(script: str, *args: str, timeout: int = 3600) -> str:
    """Run one of the R scripts beside this file; R_LIBS_USER is taken from
    the environment. Returns what it printed."""
    done = subprocess.run([RSCRIPT, os.path.join(HERE, script), *map(str, args)],
                          capture_output=True, text=True, timeout=timeout)
    if done.returncode:
        raise RuntimeError(f"{script} failed:\n{done.stdout[-2000:]}\n{done.stderr[-4000:]}")
    return done.stdout


# The experiment is unpublished and this repository is public: its groups
# are called `Case` and `Ctrl` here. `tests/real_names.json`, never committed,
# says what they are called in the file names (see tests/test_real.py).
_NAMES = os.path.join(HERE, "..", "real_names.json")
_WORDS = (__import__("json").load(open(_NAMES, encoding="utf-8"))
          if os.path.exists(_NAMES) else {})


def neutral(name: str) -> str:
    """A name from the files, as this repository writes it."""
    for word, real in _WORDS.items():
        name = name.replace(real, word)
    return name


def has_real() -> bool:
    return bool(DATA) and os.path.isdir(DATA)


def real_session():
    """The experiment as the program opens it: 3 days x Case/Ctrl x 5
    animals, gates from the FlowJo workspace. The two repeated
    acquisitions (a ` in the name) are left out. Nothing is saved."""
    from milc.session import Session
    s = Session()
    s.open_folder(DATA)
    wsp = [f for f in os.listdir(DATA) if f.endswith(".wsp")]
    if wsp:
        s.set_workspace(os.path.join(DATA, wsp[0]))
        s.set_extend_gates(True)
    s.annotate(PATTERN)
    for e in s.entries:
        if "Group" in e.meta:
            e.meta["Group"] = neutral(e.meta["Group"])
    entries = [e for e in s.entries if e.gated and not e.problem and "`" not in e.name]
    return s, entries


def real_events(population: str = CD45, channels: list[str] | None = None):
    """(biex-transformed events x channels, per event its sample index,
    the entries) of `population` in every sample of the experiment."""
    from milc.transforms import Biex
    s, entries = real_session()
    channels = channels or MARKERS
    t = Biex()
    parts, which = [], []
    for i, e in enumerate(entries):
        sm = s.load(e)
        cols = [sm.fcs.index(c) for c in channels]
        x = sm.compensated[sm.mask(population)][:, cols]
        parts.append(t(x.astype(np.float64)))
        which.append(np.full(len(x), i))
    return np.vstack(parts), np.concatenate(which), entries


def table(header: list[str], rows: list[list]) -> str:
    """A Markdown table."""
    def cell(v):
        if isinstance(v, (float, np.floating)):
            if v != v:
                return "—"
            a = abs(v)
            return f"{v:.3g}" if (a and (a < 1e-3 or a >= 1e5)) else f"{v:.4g}"
        return str(v)
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(cell(v) for v in r) + " |" for r in rows]
    return "\n".join(lines)


class Report:
    """Collects what a comparison measured: printed as it goes, and left
    as Markdown and JSON in the work folder."""

    def __init__(self, name: str):
        self.name, self.lines, self.numbers = name, [], {}

    def say(self, text: str = "") -> None:
        print(text, flush=True)
        self.lines.append(text)

    def table(self, title: str, header: list[str], rows: list[list]) -> None:
        self.say(f"\n**{title}**\n")
        self.say(table(header, rows))
        self.numbers[title] = {"header": header, "rows": rows}

    def close(self) -> None:
        folder = work(self.name)
        with open(os.path.join(folder, "report.md"), "w") as f:
            f.write("\n".join(self.lines) + "\n")
        with open(os.path.join(folder, "report.json"), "w") as f:
            json.dump(self.numbers, f, indent=1, default=float)
        print(f"\nwritten to {folder}/report.md")
