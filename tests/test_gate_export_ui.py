"""File ▸ Export gates…, driven without a screen: gates drawn in the
window are written, read back, and hold the same events."""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

QtWidgets = pytest.importorskip("PyQt6.QtWidgets")

from milc.fcs import write_fcs  # noqa: E402


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


@pytest.fixture
def folder(tmp_path):
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 600), ("3Ctrl2", 800), ("3Case1", 900), ("3Case2", 1500)):
        pop = rng.integers(0, 3, n)
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", "Foxp3"],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


def test_export_gates_from_the_window(app, folder, tmp_path):
    from milc import wsp
    from milc.gate_export import GATINGML, WSP
    from milc.gating import OR
    from milc.sample import Sample
    from milc.transforms import Biex
    from milc.ui import worker
    from milc.ui.gates_export_dialog import GatesExportDialog, GatesReportDialog
    from milc.ui.shell import MainShell

    w = MainShell()
    w.open_folder(folder)
    w.show()
    file_menu = w.menuBar().actions()[0].menu()
    assert "Export gates…" in [a.text() for a in file_menu.actions()]
    # nothing gated yet: it says so, and asks nothing
    w.choose_gates_export()
    assert w.statusBar().currentMessage() == "No sample has gates to export"

    g = w.gates
    w.tabs.setCurrentWidget(g)
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    # gates drawn as a user draws them: a rectangle, an ellipse and quadrants on
    # biexponential axes, and a combination — every sample gets them
    g.set_axes("FSC-A", "Comp-PE-A")
    g.add_rectangle(75000, 0.0, 110000, 1.0, "Cells")
    g.select("Cells")
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    t = Biex()
    cx = float(np.median(t(g.sample.column("Comp-PE-A"))))
    g.add_ellipse(cx, 0.5, 0.2, 0.15, "Round")
    corners = g.add_quadrants(cx, 0.5)
    g.add_combination(OR, [corners[0].name, corners[2].name], "Diagonal")
    s = w.session
    assert all(len(list(e.gates.walk())) == 8 for e in s.usable())

    # the dialog: the file's name follows the format
    d = GatesExportDialog(folder, 4, 4, w)
    assert d.format.currentData() == WSP and d.path.text() == os.path.join(folder, "gates.wsp")
    assert "not been opened in FlowJo" in d.about.text() and d.counts.isChecked()
    d.format.setCurrentIndex(d.format.findData(GATINGML))
    assert d.path.text().endswith("gates.xml") and "ISAC" in d.about.text()
    d.format.setCurrentIndex(d.format.findData(WSP))
    assert d.path.text().endswith("gates.wsp")
    d.close()

    # the workspace: what is read back holds, event for event, what the window shows
    out = str(tmp_path / "from_window.wsp")
    report = w.export_gates(out)
    assert report.files == [out] and report.samples == 4 and report.written == 4 * 7
    assert not report.skipped and w.statusBar().currentMessage().startswith("Wrote from_window.wsp: 28")
    back = wsp.read_wsp(out)
    assert back.skipped == [] and len(back.samples) == 4
    for e in s.usable():
        own = s.load(e).masks()
        read = back.sample(os.path.basename(e.path))
        again = Sample.open(e.path, compensation=read.compensation, gates=read.gates).masks()
        assert set(again) == set(own)
        assert all(np.array_equal(own[p], again[p]) for p in own), e.name
        assert read.gates.find("Cells/Round").reference_count == int(own["Cells/Round"].sum()) > 20
    shown = GatesReportDialog(report, w)
    assert shown.summary.text() == report.summary() and shown.text.toPlainText() == report.text()
    shown.close()

    # only the samples ticked in Samples; Gating-ML by the name of the file
    s.usable()[0].include = False
    xml = str(tmp_path / "from_window.xml")
    report = w.export_gates(xml, ticked_only=True, counts=False)
    assert report.files == [xml] and report.samples == 3 and report.format == GATINGML
    assert "Logicle nearest to it" in report.text()
    listed = [x.get("name") for x in ET.parse(xml).getroot().iter("sample")]
    assert len(listed) == 3 and os.path.basename(s.usable()[0].path) not in listed
    s.usable()[0].include = True

    # a population made of a cluster cannot be written: it is in the report, the rest is written
    from milc.gating import ClusterGate, GateNode
    for e in s.usable():
        e.gates.find("Cells").add(GateNode("Cluster 3", ClusterGate(
            ["Comp-FITC-A"], [t], [1.0], np.zeros((2, 1)), [0])))
        e.sample = None
    report = w.export_gates(str(tmp_path / "less.wsp"))
    assert report.written == 28 and len(report.skipped) == 4
    shown = GatesReportDialog(report, w)
    said = shown.text.toPlainText()
    assert "Not written:" in said and "Cells/Cluster 3 (4 samples)" in said
    assert "4 not written" in shown.summary.text()
    shown.close()
    w.close()
