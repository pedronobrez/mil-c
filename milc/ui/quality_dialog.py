"""
The quality of one acquisition, to look at before anything is removed:
the flow rate over time and, for one population of one parameter, its
median in each bin of events and its slow level, with the stretches
`qc.check` would remove shaded. The limits it checks with are here to
change, and to check again with.

Nothing is removed by looking. Applying puts a population, "Clean", at
the top of the gating tree — of this sample, or of every sample ticked in
a table that says first how much each would lose.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .help_window import describe
from .. import qc, qc_plots
from ..gating import time_node
from . import theme, worker
from .flow import captioned, flow

COLUMNS = ["Sample", "Events", "To remove", "Flow rate", "Signal", "Drift", "Shift", "Note"]
LIMITS = ("mad", "consecutive", "events_per_bin", "drift", "shift")


def limits_text(settings: dict) -> str:
    """The limits of a check, in the words of the controls."""
    bins = settings.get("events_per_bin")
    drift, shift = settings.get("drift"), settings.get("shift")
    return (f"MAD limit {settings.get('mad'):g}, {settings.get('consecutive')} consecutive "
            f"bins, {f'{bins:,d}' if bins else 'auto'} events per bin, "
            + (f"drift limit {drift:g}" if drift else "drift not judged")
            # a project applied before the short pass existed did not judge it
            + (f", shift limit {shift:g}" if shift else ", short shifts not judged"))


class QualityDialog(QtWidgets.QDialog):
    def __init__(self, tab, parent=None):
        super().__init__(parent)
        describe(self, "quality-control")
        self.tab, self.session, self.entry = tab, tab.session, tab.entry
        self.setWindowTitle(f"Quality of the acquisition: {self.entry.name}")
        self._controls()
        #: the limits `result` was found with, and those of the survey
        self.used = self.settings()
        self.survey_used = dict(self.used)
        self.result = tab.check_quality(**self.used)
        #: (entry, result) for every sample, once they have been checked
        self.survey_results: list = []
        # FCS 3 says how many seconds one unit of Time is; without it the
        # axis is in the file's units and says so
        try:
            self.step = float(tab.sample.fcs.keywords.get("$TIMESTEP", "") or 0)
        except ValueError:
            self.step = 0.0

        self.summary = QtWidgets.QLabel(self.result.summary)
        self.summary.setWordWrap(True)
        self.summary.setProperty("role", "strong" if self.result.checked else "warning")
        self.channel = QtWidgets.QComboBox()
        self._fill_channels()
        note = QtWidgets.QLabel(
            "Each point is one bin of events, in the order they were acquired. Shaded: the "
            "stretches to remove. Dashed: how far the flow rate, and the running median of "
            "the signal, may go before a bin is suspect; the bins beyond are marked. Thick "
            "line: the slow level of the signal, red over a stretch that went beyond the "
            "dotted limits: a drift, taken from where it nears them. Red bar: a short stretch "
            "whose level stood apart from its neighbours' (dashed bar): a shift. A parameter "
            "with several populations is followed one population at a time.")
        note.setWordWrap(True)
        note.setProperty("role", "caption")
        self.rate_plot, self.signal_plot = self._plot(), self._plot()
        self.signal_plot.setXLink(self.rate_plot)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setVisible(False)
        self.state = QtWidgets.QLabel("")
        self.state.setProperty("role", "caption")
        self.state.setWordWrap(True)

        self.survey_button = QtWidgets.QPushButton("Check every sample")
        self.survey_button.setToolTip("How much each sample would lose, before anything is removed")
        self.here_button = QtWidgets.QPushButton("Apply to this sample")
        self.here_button.setEnabled(self.result.checked)
        self.every_button = QtWidgets.QPushButton("Apply to every sample")
        self.every_button.setToolTip("Each sample ticked in the table loses its own bad "
                                     "stretches; one not ticked keeps every event. "
                                     "Check every sample first.")
        self.every_button.setEnabled(False)
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.remove_button.setToolTip("Take “Clean” out of every sample's gating tree: "
                                      "the gates go back to where they were")
        close = QtWidgets.QPushButton("Close")
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.survey_button)
        buttons.addWidget(self.remove_button)
        buttons.addStretch(1)
        buttons.addWidget(self.here_button)
        buttons.addWidget(self.every_button)
        buttons.addWidget(close)

        top = QtWidgets.QHBoxLayout()
        cap = QtWidgets.QLabel("Signal of")
        cap.setProperty("role", "caption")
        top.addWidget(cap)
        top.addWidget(self.channel)
        top.addStretch(1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.summary)
        layout.addWidget(flow(captioned("MAD limit", self.mad),
                              captioned("Consecutive bins", self.consecutive),
                              captioned("Events per bin", self.events_per_bin),
                              captioned("Drift limit", self.drift),
                              captioned("Shift limit", self.shift), self.again_button))
        layout.addWidget(note)
        layout.addWidget(self.rate_plot, 2)
        layout.addLayout(top)
        layout.addWidget(self.signal_plot, 3)
        layout.addWidget(self.table, 2)
        layout.addWidget(self.state)
        layout.addLayout(buttons)
        self.resize(900, 760)

        self.channel.currentIndexChanged.connect(self.draw)
        for box in (self.mad, self.consecutive, self.events_per_bin, self.drift, self.shift):
            box.valueChanged.connect(self._limits_changed)
        self.again_button.clicked.connect(self.check_again)
        self.survey_button.clicked.connect(self.survey)
        self.here_button.clicked.connect(self.apply_here)
        self.every_button.clicked.connect(self.apply_everyone)
        self.remove_button.clicked.connect(self.remove)
        close.clicked.connect(self.accept)
        self._applied()
        self.draw()

    # --- the limits -------------------------------------------------------------------

    def _controls(self) -> None:
        """The limits of `qc.check`, at what was last applied or at its
        defaults."""
        import inspect
        usual = {k: v.default for k, v in inspect.signature(qc.check).parameters.items()
                 if k in LIMITS}
        start = {**usual, **{k: v for k, v in (self.session.quality_settings or {}).items()
                             if k in LIMITS}}
        self.mad = QtWidgets.QDoubleSpinBox()
        self.mad.setRange(1.0, 50.0)
        self.mad.setSingleStep(0.5)
        self.mad.setDecimals(1)
        self.mad.setValue(float(start["mad"]))
        self.mad.setToolTip("How many MADs from the acquisition's own level the signal has to "
                            "lie for a bin to be suspect, and the flow rate from the usual "
                            f"one. Larger removes less. Usually {usual['mad']:g}.")
        self.consecutive = QtWidgets.QSpinBox()
        self.consecutive.setRange(1, 100)
        self.consecutive.setValue(int(start["consecutive"]))
        self.consecutive.setToolTip("Good bins between two removed stretches are kept only "
                                    "when there are at least this many. Usually "
                                    f"{usual['consecutive']}.")
        self.events_per_bin = QtWidgets.QSpinBox()
        self.events_per_bin.setRange(0, 100000)
        self.events_per_bin.setSingleStep(50)
        self.events_per_bin.setSpecialValueText("auto")
        self.events_per_bin.setValue(int(start["events_per_bin"] or 0))
        self.events_per_bin.setToolTip(f"auto: 1/{qc.MAX_BINS} of the file and at least "
                                       f"{qc.MIN_BIN} events. Larger bins are steadier and "
                                       "see less of a short defect.")
        self.drift = QtWidgets.QDoubleSpinBox()
        self.drift.setRange(0.0, 5.0)
        self.drift.setSingleStep(0.05)
        self.drift.setDecimals(2)
        self.drift.setSpecialValueText("off")
        self.drift.setValue(float(start["drift"] or 0))
        self.drift.setToolTip("The smallest lasting shift of a population's level that is "
                              "taken for a drift, as a share of the spread of the "
                              "parameter's events. Smaller finds smaller drifts; off does "
                              f"not look for them. Usually {usual['drift']:g}.")
        self.shift = QtWidgets.QDoubleSpinBox()
        self.shift.setRange(0.0, 5.0)
        self.shift.setSingleStep(0.05)
        self.shift.setDecimals(2)
        self.shift.setSpecialValueText("off")
        self.shift.setValue(float(start["shift"] or 0))
        self.shift.setToolTip("The smallest move of a short stretch (a few bins, up to the "
                              "drift's window) away from the bins on either side that is "
                              "taken for a shift, as a share of the spread of the "
                              "parameter's events; it has to stand out by the MAD limit "
                              "too. Smaller finds smaller shifts; off does not look for "
                              f"them. Usually {usual['shift']:g}.")
        for box, wide in ((self.mad, 76), (self.consecutive, 64), (self.events_per_bin, 96),
                          (self.drift, 88), (self.shift, 88)):
            box.setFixedWidth(wide)
        self.again_button = QtWidgets.QPushButton("Check again")
        self.again_button.setToolTip("Check this sample with these limits — and every "
                                     "sample, where they have been checked")
        self.again_button.setEnabled(False)

    def settings(self) -> dict:
        """The limits as the controls have them, to call `qc.check` with."""
        return {"mad": float(self.mad.value()), "consecutive": int(self.consecutive.value()),
                "events_per_bin": int(self.events_per_bin.value()) or None,
                "drift": round(float(self.drift.value()), 2),
                "shift": round(float(self.shift.value()), 2)}

    def _limits_changed(self) -> None:
        self.again_button.setEnabled(self.settings() != self.used)

    def check_again(self) -> None:
        """This sample with the limits as they are now; and every sample,
        where the table is of other limits."""
        self.used = self.settings()
        self.result = self.tab.check_quality(**self.used)
        self.summary.setText(self.result.summary)
        self.summary.setProperty("role", "strong" if self.result.checked else "warning")
        self.summary.style().polish(self.summary)
        self.here_button.setEnabled(self.result.checked)
        self.again_button.setEnabled(False)
        self._fill_channels()
        self.draw()
        if self.survey_results and self.survey_used != self.used:
            self.survey()

    def _fill_channels(self) -> None:
        """The populations followed, on the one that was being looked at —
        or, at first, on the one that went wrong for longest."""
        was = self.channel.currentText()
        self.channel.blockSignals(True)
        self.channel.clear()
        self.channel.addItems(self.result.populations)
        if was in self.result.populations:
            self.channel.setCurrentIndex(self.result.populations.index(was))
        elif self.result.checked:
            b = self.result.bins
            bad = b.bad_signal | b.bad_drift | b.bad_shift
            if bad.any():
                self.channel.setCurrentIndex(int(bad.sum(axis=0).argmax()))
        self.channel.blockSignals(False)

    @staticmethod
    def _plot() -> pg.PlotWidget:
        plot = pg.PlotWidget()
        plot.setBackground(theme.background())
        plot.setMenuEnabled(False)
        plot.hideButtons()
        plot.setMouseEnabled(True, False)
        return plot

    # --- the plots ------------------------------------------------------------------

    def set_channel(self, name: str) -> None:
        """Look at this population, named as `qc.Result.populations` names it."""
        self.channel.setCurrentIndex(self.result.populations.index(name))

    def series(self) -> dict:
        """What the two plots are drawn from, in the units of their axes:
        the time of each bin, its flow rate, the median of the population
        chosen (not a number where a bin has none of it), the limits of
        both, the slow level with its own limits, the short stretches that
        stood apart and the stretches shaded."""
        return qc_plots.series(self.result, self.channel.currentIndex(), self.step)

    def draw(self) -> None:
        for plot in (self.rate_plot, self.signal_plot):
            plot.clear()
            theme.style_axes(plot.plotItem)
        if not self.result.checked:
            return
        s = self.series()
        unit = "s" if self.step else "as the file stores it"
        self.rate_plot.setLabel("left", "Events per second" if self.step
                                else "Events per unit of time")
        self.signal_plot.setLabel("left", self.channel.currentText())
        self.signal_plot.setLabel("bottom", f"Time ({unit})")
        for plot, y, bad, limits in (
                (self.rate_plot, s["rate"], s["bad_rate"], s["rate_limits"]),
                (self.signal_plot, s["median"], s["bad_signal"], s["limits"])):
            shade = QtGui.QColor(theme.danger())
            shade.setAlpha(46)
            for lo, hi in s["shaded"]:
                region = pg.LinearRegionItem((lo, hi), movable=False, brush=pg.mkBrush(shade),
                                             pen=pg.mkPen(None))
                region.setZValue(-10)
                plot.addItem(region)
            for level in limits:
                plot.addItem(pg.InfiniteLine(level, angle=0, pen=pg.mkPen(
                    theme.muted(), width=1, style=QtCore.Qt.PenStyle.DashLine)))
            plot.addItem(pg.PlotCurveItem(s["time"], y, connect="finite",
                                          pen=pg.mkPen(theme.accent(), width=1.2)))
            bad = bad & np.isfinite(y)
            plot.addItem(pg.ScatterPlotItem(s["time"][bad], y[bad], size=6, pen=None,
                                            brush=pg.mkBrush(theme.danger())))
            if plot is self.signal_plot:
                for level in s["slow_limits"]:
                    plot.addItem(pg.InfiniteLine(level, angle=0, pen=pg.mkPen(
                        theme.muted(), width=1, style=QtCore.Qt.PenStyle.DotLine)))
                plot.addItem(pg.PlotCurveItem(s["time"], s["slow"], connect="finite",
                                              pen=pg.mkPen(theme.foreground(), width=2.2)))
                plot.addItem(pg.PlotCurveItem(
                    s["time"], np.where(s["bad_drift"], s["slow"], np.nan), connect="finite",
                    pen=pg.mkPen(theme.danger(), width=3.2)))
                for lo, hi, level, near in s["shifts"]:
                    plot.addItem(pg.PlotCurveItem([lo, hi], [near, near], pen=pg.mkPen(
                        theme.foreground(), width=1.2, style=QtCore.Qt.PenStyle.DashLine)))
                    plot.addItem(pg.PlotCurveItem([lo, hi], [level, level],
                                                  pen=pg.mkPen(theme.danger(), width=3.2)))
                limits = [min(limits[0], s["slow_limits"][0]), max(limits[1], s["slow_limits"][1])]
            lo, hi = min(float(np.nanmin(y)), limits[0]), max(float(np.nanmax(y)), limits[1])
            plot.setYRange(lo, hi, padding=0.06)
        self.rate_plot.setXRange(float(s["time"][0]), float(s["time"][-1]), padding=0.02)

    # --- every sample ----------------------------------------------------------------

    def survey(self) -> None:
        """Check every sample off the window's thread and list what each
        would lose."""
        entries = self.session.usable()
        used = dict(self.used)
        self.survey_results, self.survey_used = [], used
        self.survey_button.setEnabled(False)
        self.every_button.setEnabled(False)
        self.table.setVisible(True)
        self.table.setRowCount(0)
        self.state.setText("Checking every sample…")

        def work(progress):
            out = []
            for i, e in enumerate(entries):
                out.append((e, self.session.quality(e, **used)))
                progress(i + 1, len(entries), e.name)
            return out

        worker.run(work, on_done=self._surveyed, on_failed=self._failed,
                   on_progress=lambda i, n, name: self.state.setText(f"Checked {name} ({i} of {n})"))

    def _failed(self, message: str) -> None:
        self.state.setText(message)
        self.survey_button.setEnabled(True)

    def _surveyed(self, results: list) -> None:
        self.survey_results = results
        self.table.setRowCount(len(results))
        right = QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter
        for row, (e, r) in enumerate(results):
            name = QtWidgets.QTableWidgetItem(e.name)
            if r.checked:
                name.setFlags(name.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                name.setCheckState(QtCore.Qt.CheckState.Checked)
            cells = [f"{e.events:,d}",
                     f"{100 * r.removed:.1f}%" if r.checked else "—",
                     f"{r.counts[qc.RATE]:,d}" if r.checked else "—",
                     f"{r.counts[qc.SIGNAL]:,d}" if r.checked else "—",
                     f"{r.counts[qc.DRIFT]:,d}" if r.checked else "—",
                     f"{r.counts[qc.SHIFT]:,d}" if r.checked else "—",
                     ", ".join(r.blamed) if r.checked else r.reason]
            self.table.setItem(row, 0, name)
            for col, text in enumerate(cells, start=1):
                item = QtWidgets.QTableWidgetItem(text)
                if col < len(cells):
                    item.setTextAlignment(right)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
        self.survey_button.setEnabled(True)
        self.every_button.setEnabled(any(r.checked for _, r in results))
        lost = sum(int((~r.keep).sum()) for _, r in results)
        total = sum(len(r.keep) for _, r in results)
        self.state.setText(f"{len(results)} samples would lose {lost:,d} of {total:,d} events "
                           "in all. Untick a sample to leave it as it is.")

    def ticked(self) -> list[str]:
        return [self.table.item(row, 0).text() for row in range(self.table.rowCount())
                if self.table.item(row, 0).checkState() == QtCore.Qt.CheckState.Checked]

    def set_ticked(self, name: str, on: bool) -> None:
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item.text() == name and item.flags() & QtCore.Qt.ItemFlag.ItemIsUserCheckable:
                item.setCheckState(QtCore.Qt.CheckState.Checked if on
                                   else QtCore.Qt.CheckState.Unchecked)

    # --- applying ----------------------------------------------------------------------

    def apply_here(self) -> None:
        if self.result.checked:
            self.tab.apply_quality([(self.entry, self.result.intervals)], self.used)
            self._applied()

    def apply_everyone(self) -> None:
        """Every sample ticked loses its own stretches; one not ticked
        keeps every event."""
        ticked = set(self.ticked())
        found = [(e, r.intervals if e.name in ticked else [])
                 for e, r in self.survey_results if r.checked]
        if found:
            self.tab.apply_quality(found, self.survey_used)
            self._applied()

    def remove(self) -> None:
        self.tab.remove_quality()
        self._applied()

    def _applied(self) -> None:
        """Say what the trees hold now."""
        cleaned = [(e, time_node(e.gates)) for e in self.session.usable()]
        cleaned = [(e, node) for e, node in cleaned if node is not None]
        self.remove_button.setEnabled(bool(cleaned))
        if not cleaned:
            self.state.setText("Nothing is removed from any sample.")
            return
        losing = sum(bool(node.gate.intervals) for _, node in cleaned)
        name = cleaned[0][1].name
        self.state.setText(
            f"“{name}” is the first population of every sample, and every gate is under "
            f"it. {losing} of {len(cleaned)} samples have stretches of time removed. "
            + (f"Last applied with {limits_text(self.session.quality_settings)}."
               if self.session.quality_settings else ""))
