"""
Against the acquisitions in Data/ and the FlowJo workspace beside them.

Measured on 2026-11-22, 12 samples, 72 gates: FlowJo evaluates polygon
gates on a 256-bin grid (`gateResolution`) by a rule that was not
reproduced — four guesses were tried and none matched — so ours are exact
geometry and differ from its counts by at most 0.77% on any gate, and by
at most 3 events on the final population. The tolerances below are that
measurement with a margin, not a target.
"""

import json
import os

import pytest

from milc.sample import Sample
from milc.wsp import read_wsp

DATA = os.path.join(os.path.dirname(__file__), "..", "Data")

# The experiment is unpublished, and this repository is public: its samples
# are called by neutral names here (`3Case2`, `14Ctrl4`; the FlowJo files
# `workspace`). `real_names.json` beside this file, never committed, says
# what each neutral word is in the names of the files:
#     {"Case": "...", "Ctrl": "...", "workspace": "..."}
_NAMES = os.path.join(os.path.dirname(__file__), "real_names.json")
_WORDS = json.load(open(_NAMES, encoding="utf-8")) if os.path.exists(_NAMES) else {}


def real(name: str) -> str:
    """A name as written here, as the files on disk have it."""
    for neutral, word in _WORDS.items():
        name = name.replace(neutral, word)
    return name


WSP = os.path.join(DATA, real("workspace") + ".wsp")
TCRB = "Cells/Single Cells/Single Cells/Live/CD45/Tcrb"

pytestmark = [pytest.mark.real,
              pytest.mark.skipif(not os.path.exists(WSP), reason="no real data here")]


@pytest.fixture(scope="module")
def workspace():
    return read_wsp(WSP)


def test_workspace_structure(workspace):
    gated = [s for s in workspace.samples if s.compensation is not None]
    assert len(gated) == 12
    s = workspace.sample(real("3Case2.fcs"))
    assert [n.path for n in s.gates.leaves()] == [TCRB]
    assert s.gates.find(TCRB).reference_count == 638
    # the matrix is the one edited in FlowJo, not the file's own
    i, j = s.compensation.channels.index("FITC-A"), s.compensation.channels.index("PE-A")
    assert s.compensation.matrix[i, j] == pytest.approx(0.15997)


def test_counts_against_flowjo(workspace):
    for w in workspace.samples:
        if w.compensation is None:
            continue
        s = Sample.open(w.locate(DATA), compensation=w.compensation, gates=w.gates)
        assert len(s) == w.gates.reference_count
        for node in w.gates.walk():
            if node.gate is None:
                continue
            ours, ref = s.count(node.path), node.reference_count
            assert abs(ours - ref) <= max(0.012 * ref, 3), (w.name, node.path)
        assert abs(s.count(TCRB) - w.gates.find(TCRB).reference_count) <= 5


# --- the analysis, against what FlowJo's plugins wrote -----------------------
#
# Measured on 2026-11-22 on FlowJo's two concatenated files (2040 and 2100
# events, five markers):
#   - the medians of Biex() x 1023 per FlowJo metacluster are FlowJo's own
#     PopMFI table to within 0.69 of a channel: that is the transform, and
#     that is what the plugins are fed;
#   - our FlowSOM against FlowJo's labels: adjusted Rand 0.74-0.91 over 8
#     seeds (mean 0.86 and 0.81). Our own runs against each other: mean 0.81
#     and 0.79 — so FlowJo's run is as close to ours as ours are to each other;
#   - neighbours kept (k=15): UMAP 0.470 / 0.479 ours, 0.472 / 0.477 FlowJo;
#     tSNE 0.588 ours, 0.580 FlowJo.

import csv

import numpy as np

from milc.dataset import Dataset
from milc.transforms import Biex

RUNS = os.path.join(DATA, real("workspace"))
FIVE = ["FITC-A", "PE-A", "APC-A", "BV421-A", "BV510-A"]
CASES = [("concat_pedro_2_1.fcs", "EF93", "OSQG"),
         ("concat_UMAP_group_1.fcs", "G7WJ", "IQGK")]


def _flowjo(name, som, um):
    ds = Dataset.from_fcs(os.path.join(DATA, name), sample_column=None)
    X = ds.matrix(FIVE, Biex(), 1023)
    d = os.path.join(RUNS, "FlowSOM", f"Ungated_{som}")
    labels = np.loadtxt(os.path.join(d, f"FlowSOM_Results_{som}.csv"),
                        delimiter=",", skiprows=1)[:, 1].astype(int)
    with open(os.path.join(d, f"MetaCluster_MFI_{som}_.PopMFI.csv")) as f:
        mfi = np.array([[float(v) for v in r[1:]] for r in list(csv.reader(f))[1:]])
    U = np.loadtxt(os.path.join(RUNS, "UMAP", f"UMAP_RunID_{um}", f"UMAP_Results_{um}.csv"),
                   delimiter=",", skiprows=1)[:, 1:]
    return X, labels, mfi, U


@pytest.mark.parametrize("case", CASES)
def test_biex_is_what_flowjo_feeds_its_plugins(case):
    X, labels, mfi, _ = _flowjo(*case)
    ours = np.array([np.median(X[labels == i], axis=0) for i in range(1, 6)])
    assert np.abs(ours - mfi).max() < 1.0


@pytest.mark.parametrize("case", CASES)
def test_flowsom_agrees_with_flowjo(case):
    from milc import flowsom as fs
    X, labels, _, _ = _flowjo(*case)
    ari = [fs.adjusted_rand(fs.flowsom(X, 5, 10, 10, seed=s).labels, labels)
           for s in range(4)]
    assert np.mean(ari) > 0.7 and min(ari) > 0.6


@pytest.mark.parametrize("case", CASES)
def test_umap_keeps_neighbours_as_flowjos_does(case):
    from milc import embed
    X, _, _, U = _flowjo(*case)
    theirs = embed.neighbour_preservation(X, U)
    ours = embed.neighbour_preservation(X, embed.umap(X, seed=0))
    assert ours > theirs - 0.03


def test_files_outside_the_workspace_get_the_shared_gates(workspace):
    """Measured 2026-11-22: the 12 samples of the workspace have one matrix
    and one gating tree between them, identical to the group's. The 20
    files it does not list are given both."""
    from milc.session import Session
    gates, matrix = workspace.shared()
    assert matrix.name == "Acquisition-defined-copy"
    assert [n.path for n in gates.leaves()] == [TCRB]
    assert all(n.reference_count is None for n in gates.walk())

    s = Session()
    s.open_folder(DATA)
    assert s.can_extend() and s.population == TCRB
    usable = [e for e in s.entries if s.available(e)]
    assert len(usable) == 32 and sum(e.borrowed for e in usable) == 20
    by_name = {e.name: e for e in s.entries}
    assert s.count(by_name[real("3Ctrl1")]) == 109 and s.count(by_name[real("3Ctrl2")]) == 175
    # whatever the project in Data/ has ticked: here, in memory, every file
    for e in s.entries:
        e.include = True
    assert sorted(map(sorted, s.repeats())) == sorted(
        [sorted(map(real, pair)) for pair in (["14Case4", "14Case4`"], ["7Case3", "7Case3`"])])

    s.set_extend_gates(False)
    assert sum(s.available(e) for e in s.entries) == 12
