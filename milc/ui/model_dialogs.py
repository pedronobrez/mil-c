"""
The three models of the Gates tab's "Model…" button — proliferation from
a dye's dilution, the cell cycle from DNA content, kinetics against Time
— each a dialog on the population shown, for every sample selected:
the data with what was fitted, the settings of the model, a table with a
row for each sample, and the ways out (a figure, a CSV, a comparison
between groups).

The models themselves are in `milc.proliferation`, `milc.cellcycle` and
`milc.kinetics`, and know nothing of windows.
"""

from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import cellcycle, compare, kinetics, proliferation
from ..stats import StatsTable
from . import palette, theme
from .settings import settings
from .statistics_tab import _Cell, _format

PROLIFERATION, CELL_CYCLE, KINETICS = "Proliferation", "Cell cycle", "Kinetics"
SETTING_FOLDER = "statistics/folder"


class Source:
    """What the Gates tab shows, as a model asks for it: the samples
    selected, the population, the parameters, and the events' values as
    the files store them (compensated, on no display scale)."""

    def __init__(self, tab):
        self.session = tab.session
        self.entries = list(tab.overlaid())
        self.entry = tab.entry
        self.path = tab.node.path
        self.population = tab.node.name if tab.node.gate is not None else "All events"
        self.parameters = [(tab.x_axis.itemText(i), tab.x_axis.itemData(i))
                           for i in range(tab.x_axis.count())]
        self.dim = tab.x_axis.currentData()
        # the column the Gates tab colours by; None: a colour for each sample
        by = tab.colour_by.currentData()
        self.colour_by = by if by in tab.session.columns else None

    def label(self, dim: str) -> str:
        return next((text for text, d in self.parameters if d == dim), dim)

    def group(self, entry) -> str:
        return entry.name if self.colour_by is None else (entry.meta.get(self.colour_by, "") or "—")

    def sample(self, entry):
        return self.session.load(entry)

    def values(self, entry, dim: str) -> np.ndarray:
        """`dim` of the population's events in this sample; ValueError
        where the sample has no such population or parameter."""
        sample = self.sample(entry)
        masks = sample.masks()
        if self.path not in masks:
            raise ValueError(f"{self.population} is not gated in this sample")
        try:
            return sample.column(dim)[masks[self.path]].astype(np.float64)
        except (ValueError, KeyError, IndexError):
            raise ValueError(f"no parameter {dim} in this sample") from None

    def time(self, entry) -> tuple[np.ndarray, float]:
        """(Time of the population's events, seconds in one of its units
        or 0 where the file does not say)."""
        sample = self.sample(entry)
        channel = kinetics.time_channel(sample.fcs.channels)
        if channel is None:
            raise ValueError("the file has no Time parameter")
        return self.values(entry, channel), kinetics.timestep(sample.fcs.keywords)


def _spin(value: float, lo: float, hi: float, decimals: int = 2, step: float = 0.1,
          tip: str = "") -> QtWidgets.QDoubleSpinBox:
    box = QtWidgets.QDoubleSpinBox()
    box.setDecimals(decimals)
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setValue(value)
    box.setKeyboardTracking(False)
    box.setToolTip(tip)
    return box


class _ModelDialog(QtWidgets.QDialog):
    """The plot, the settings, the table and the ways out. A model fills
    in `_controls`, `columns`, `_fit` and `_draw`."""

    KIND = ""
    NAME = "model"
    #: the column a comparison between groups starts on
    COMPARE = ""
    DEFINITIONS = ""

    def __init__(self, source: Source, parent=None):
        super().__init__(parent)
        self.source = source
        self.setWindowTitle(f"{self.KIND} — {source.population}")
        #: sample → what was fitted, and → why it could not be
        self.fits: dict = {}
        self.problems: dict = {}
        self.result: StatsTable | None = None
        self.shown: str | None = None
        self._filling = False
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self.compute)

        self.parameter = QtWidgets.QComboBox()
        for text, dim in source.parameters:
            self.parameter.addItem(text, dim)
        self.form = QtWidgets.QFormLayout()
        self.form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self.form.addRow("Parameter", self.parameter)
        self._controls()
        self.parameter.currentIndexChanged.connect(self.changed)
        panel = QtWidgets.QWidget()
        side = QtWidgets.QVBoxLayout(panel)
        side.setContentsMargins(0, 0, 0, 0)
        side.addLayout(self.form)
        side.addStretch(1)
        panel.setFixedWidth(330)

        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.setMinimumHeight(260)
        theme.style_axes(self.plot)
        self.said = QtWidgets.QLabel("")
        self.said.setProperty("role", "strong")
        self.said.setWordWrap(True)
        right = QtWidgets.QVBoxLayout()
        right.addWidget(self.said)
        right.addWidget(self.plot, 1)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(panel)
        top.addLayout(right, 1)

        self.table = QtWidgets.QTableWidget(0, 0)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setMinimumHeight(150)
        self.table.currentCellChanged.connect(self._row_chosen)
        self.caption = QtWidgets.QLabel(self.DEFINITIONS)
        self.caption.setProperty("role", "caption")
        self.caption.setWordWrap(True)

        self.figure_button = QtWidgets.QPushButton("Figure…")
        self.compare_button = QtWidgets.QPushButton("Compare groups…")
        self.compare_button.setToolTip("A test of one column of the table between groups of "
                                       "samples: each sample is one observation")
        self.export_button = QtWidgets.QPushButton("Export CSV…")
        self.copy_button = QtWidgets.QPushButton("Copy")
        self.copy_button.setToolTip("The table, to paste into Excel or Prism")
        self.close_button = QtWidgets.QPushButton("Close")
        self.figure_button.clicked.connect(self.make_figure)
        self.compare_button.clicked.connect(self._ask_comparison)
        self.export_button.clicked.connect(self._ask_export)
        self.copy_button.clicked.connect(self.copy)
        self.close_button.clicked.connect(self.accept)
        buttons = QtWidgets.QHBoxLayout()
        for b in (self.figure_button, self.compare_button, self.export_button, self.copy_button):
            buttons.addWidget(b)
        buttons.addStretch(1)
        buttons.addWidget(self.close_button)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(top, 3)
        layout.addWidget(self.table, 2)
        layout.addWidget(self.caption)
        layout.addLayout(buttons)
        self.resize(1120, 780)

    # --- what a model says of itself ------------------------------------------------

    def _controls(self) -> None:
        raise NotImplementedError

    def columns(self) -> list[str]:
        raise NotImplementedError

    def _fit(self, entry):
        """(what was fitted, the numbers of its row)."""
        raise NotImplementedError

    def _prepare(self) -> None:
        """Before the samples are fitted: what they share."""

    def _draw(self) -> None:
        raise NotImplementedError

    def content(self):
        """What is on the plot, for a figure; None when there is nothing."""
        raise NotImplementedError

    # --- fitting ------------------------------------------------------------------------

    def changed(self, *_) -> None:
        if not self._filling:
            self._timer.start()

    def set_parameter(self, dim: str) -> None:
        i = self.parameter.findData(dim)
        if i >= 0:
            self.parameter.setCurrentIndex(i)
        self.compute()

    def compute(self) -> None:
        """Fit every sample selected, with the settings as they are."""
        self._timer.stop()
        s = self.source
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            self.fits, self.problems = {}, {}
            names = self.columns()
            table = StatsTable(["Sample", *s.session.columns, "Population", "Path", *names,
                                "Note"])
            try:
                self._prepare()
                shared = ""
            except ValueError as exc:
                shared = str(exc)
            for e in s.entries:
                row = [float("nan")] * len(names)
                try:
                    if shared:
                        raise ValueError(shared)
                    self.fits[e.name], row = self._fit(e)
                except ValueError as exc:
                    self.problems[e.name] = str(exc)
                table.rows.append([e.name, *[e.meta.get(c, "") for c in s.session.columns],
                                   s.population, s.path, *row, self.problems.get(e.name, "")])
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._show_table(table)

    def _show_table(self, table: StatsTable) -> None:
        self.result = table
        hidden = {"Path"} | ({"Note"} if not self.problems else set())
        shown = [i for i, c in enumerate(table.columns) if c not in hidden]
        self._filling = True
        self.table.setSortingEnabled(False)
        self.table.clear()
        self.table.setColumnCount(len(shown))
        self.table.setHorizontalHeaderLabels([table.columns[i] for i in shown])
        self.table.setRowCount(len(table.rows))
        for r, row in enumerate(table.rows):
            for k, i in enumerate(shown):
                self.table.setItem(r, k, _Cell(row[i], _format(table.columns[i], row[i])))
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        self._filling = False
        names = [e.name for e in self.source.entries]
        keep = self.shown if self.shown in names else (
            self.source.entry.name if self.source.entry is not None
            and self.source.entry.name in names else (names[0] if names else None))
        self.show_sample(keep)
        usable = bool(self.fits)
        for b in (self.figure_button, self.compare_button):
            b.setEnabled(usable)

    def _row_chosen(self, row: int, *_) -> None:
        if self._filling or row < 0 or self.table.item(row, 0) is None:
            return
        self.show_sample(self.table.item(row, 0).text(), select=False)

    def show_sample(self, name: str | None, select: bool = True) -> None:
        """Put this sample's fit on the plot."""
        self.shown = name
        if select and name is not None:
            self._filling = True
            for r in range(self.table.rowCount()):
                if self.table.item(r, 0).text() == name:
                    self.table.setCurrentCell(r, 0)
            self._filling = False
        self.plot.clear()
        if self.plot.plotItem.legend is not None:
            self.plot.plotItem.legend.clear()
        self._draw()

    def _legend(self):
        legend = self.plot.plotItem.legend or self.plot.plotItem.addLegend(offset=(-8, 8))
        legend.clear()
        legend.setBrush(theme.legend_brush())
        legend.setPen(theme.legend_pen())
        legend.setLabelTextColor(theme.foreground())
        return legend

    def _draw_fit(self, content, linear_ticks: bool) -> None:
        """A histogram with the components of its model."""
        x, width = content.x, float(content.x[1] - content.x[0])
        edges = np.concatenate([x - width / 2, x[-1:] + width / 2])
        grey = QtGui.QColor(theme.muted())
        grey.setAlpha(70)
        self.plot.addItem(pg.PlotCurveItem(edges, content.counts, stepMode="center",
                                           fillLevel=0.0, brush=pg.mkBrush(grey),
                                           pen=pg.mkPen(theme.muted(), width=1)))
        legend = self._legend()
        colours = palette.categorical(len(content.components))
        for comp, colour in zip(content.components, colours):
            fill = QtGui.QColor(colour)
            fill.setAlpha(70)
            item = pg.PlotCurveItem(x, comp.y, fillLevel=0.0, brush=pg.mkBrush(fill),
                                    pen=pg.mkPen(colour, width=1.4))
            self.plot.addItem(item)
            legend.addItem(item, f"{comp.name} · {comp.share:.1f}%")
        model = pg.PlotCurveItem(x, content.model, pen=pg.mkPen(theme.foreground(), width=1.6))
        self.plot.addItem(model)
        legend.addItem(model, "Model")
        for at in content.fitted or ():
            self.plot.addItem(pg.InfiniteLine(at, pen=pg.mkPen(
                theme.muted(), width=1, style=QtCore.Qt.PenStyle.DotLine)))
        self.plot.getAxis("bottom").setTicks(None if linear_ticks else [content.xticks])
        self.plot.setLabel("bottom", content.xlabel)
        self.plot.setLabel("left", "Events")
        self.plot.setXRange(edges[0], edges[-1], padding=0)
        self.plot.setYRange(0, max(float(content.counts.max()), float(content.model.max())) * 1.06,
                            padding=0)

    def _nothing(self) -> None:
        why = self.problems.get(self.shown or "", "") if self.shown else "No sample selected"
        self.said.setText(f"{self.shown}: {why}" if self.shown and why else why)

    # --- the ways out ---------------------------------------------------------------------

    def make_figure(self):
        from .figure_dialog import FigureDialog
        content = self.content()
        if content is None:
            return None
        dialog = FigureDialog(content, self, self.NAME)
        dialog.exec()
        return dialog

    def value_columns(self) -> list[str]:
        """The columns of the table that hold a number for each sample."""
        cols = self.result.columns if self.result is not None else []
        return [c for c in cols[cols.index("Path") + 1:] if c != "Note"] if cols else []

    def comparison(self, value: str, group: list[str], within: str | None = None,
                   kind: str = compare.PARAMETRIC, correction: str = compare.BH,
                   pairwise: bool = False, paired_by: str | None = None):
        """(table of results, the comparisons) of one column between the
        groups of samples `group` makes, as the Statistics tab's."""
        if paired_by:
            values, names, subjects = compare.from_table(self.result, value, group, within,
                                                         paired_by)
            problem = compare.pairing_problem(values, subjects)
            if problem:
                raise ValueError(f"The samples cannot be paired by {paired_by}: {problem}.")
            found = compare.compare(values, kind, correction, pairwise, subjects)
            return compare.to_table(found, names, correction, paired_by), found
        values, names = compare.from_table(self.result, value, group, within)
        found = compare.compare(values, kind, correction, pairwise)
        return compare.to_table(found, names, correction), found

    def _ask_comparison(self) -> None:
        if self.result is None:
            return
        from .results_dialog import NOTE, ResultsDialog
        from .statistics_tab import CompareDialog
        session = self.source.session
        dialog = CompareDialog(self.value_columns(), session.columns, session.group_by, self)
        dialog.value.setCurrentIndex(max(dialog.value.findData(self.COMPARE), 0))
        if not dialog.exec() or not dialog.groups():
            return
        value = dialog.value.currentData()
        try:
            table, found = self.comparison(
                value, dialog.groups(), dialog.within.currentData() or None,
                dialog.test.currentData(), dialog.correction.currentData(),
                dialog.pairwise.isChecked(), dialog.paired_by.currentData() or None)
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Compare groups", str(exc))
            return
        note = NOTE.format(tests=sum(1 for c in found if c.p == c.p),
                           correction=dialog.correction.currentData())
        ResultsDialog(table, f"{value} between groups",
                      f"{value} of {self.source.population}, compared between "
                      f"{' + '.join(dialog.groups())}. " + note, self,
                      f"{self.NAME}_comparison").exec()

    def export(self, path: str) -> None:
        self.result.to_csv(path)

    def _ask_export(self) -> None:
        if self.result is None:
            return
        folder = settings().value(SETTING_FOLDER, "") or self.source.session.folder
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export the table", os.path.join(folder or "", self.NAME + ".csv"),
            "Comma-separated values (*.csv)")
        if path:
            self.export(path)
            settings().setValue(SETTING_FOLDER, os.path.dirname(path))

    def copy(self) -> None:
        if self.result is not None:
            QtWidgets.QApplication.clipboard().setText(self.result.to_text())


# --- proliferation --------------------------------------------------------------------------

class ProliferationDialog(_ModelDialog):
    KIND, NAME, COMPARE = PROLIFERATION, "proliferation", "Division Index"
    DEFINITIONS = (
        "FlowJo's definitions. % Divided: the share of the cells at the start that divided "
        "(precursor frequency). Division Index: divisions per cell at the start. "
        "Proliferation Index: divisions per cell that divided. Expansion Index: how many "
        "times the culture grew. Replication Index: how many times the cells that divided "
        "grew. Generations: the highest holding 0.5% of the modelled events. SE: standard "
        "error from counting — how far the same cells acquired again would move the "
        "number, not whether the peaks are Gaussian; RMS (events per bin, and % of the "
        "tallest bin) says that. Dead cells are not seen.")

    def _controls(self) -> None:
        s = self.source
        self.control = QtWidgets.QComboBox()
        self.control.addItem("Fitted in each sample", "")
        for e in s.session.usable():
            self.control.addItem(e.name, e.name)
        self.control.setToolTip(
            "A sample that was not stimulated: its peak says which peak of the others is "
            "the undivided one, which is then looked for within a third of a division of "
            "it. Without one, the brightest peak that holds cells is taken to be undivided "
            "— wrong where every cell has divided.")
        self.generations = QtWidgets.QSpinBox()
        self.generations.setRange(1, 12)
        self.generations.setValue(proliferation.MAX_GENERATIONS)
        self.generations.setToolTip("Peaks in the model after the undivided one")
        self.grow = QtWidgets.QCheckBox("Peaks widen with each generation")
        self.grow.setToolTip("Unticked, every peak has one width on the logarithmic axis. "
                             "Ticked, the width may grow by the same fraction, up to 30%, "
                             "with each generation.")
        self.background = _spin(0.0, 0.0, 1e7, 0, 50.0,
                                "The signal of cells without the dye, in the units of the "
                                "file: only the dye above it halves. 0: none.")
        self.form.addRow("Undivided peak", self.control)
        self.form.addRow("Generations, at most", self.generations)
        self.form.addRow("Autofluorescence", self.background)
        self.form.addRow(self.grow)
        self.control.currentIndexChanged.connect(self.changed)
        self.generations.valueChanged.connect(self.changed)
        self.background.valueChanged.connect(self.changed)
        self.grow.toggled.connect(self.changed)
        self.parameter.setCurrentIndex(max(self.parameter.findData(s.dim), 0))

    def columns(self) -> list[str]:
        return proliferation.columns(self.generations.value())

    def set_control(self, name: str | None) -> None:
        self.control.setCurrentIndex(max(self.control.findData(name or ""), 0))
        self.compute()

    def _prepare(self) -> None:
        self.peak = None
        name = self.control.currentData()
        if name:
            entry = next(e for e in self.source.session.usable() if e.name == name)
            try:
                self.peak = proliferation.undivided(
                    self.source.values(entry, self.parameter.currentData()))
            except ValueError as exc:
                raise ValueError(f"control {name}: {exc}") from None

    def _fit(self, entry):
        f = proliferation.fit(self.source.values(entry, self.parameter.currentData()),
                              self.peak, self.generations.value(), self.grow.isChecked(),
                              self.background.value())
        return f, proliferation.row(f, self.generations.value())

    def content(self):
        f = self.fits.get(self.shown)
        if f is None:
            return None
        return proliferation.content(f, self.source.label(self.parameter.currentData()),
                                     self.shown)

    def _draw(self) -> None:
        content = self.content()
        if content is None:
            self._nothing()
            return
        self._draw_fit(content, linear_ticks=False)
        f = self.fits[self.shown]
        if self.peak is not None:
            self.plot.addItem(pg.InfiniteLine(self.peak.position, pen=pg.mkPen(
                theme.accent(), width=1, style=QtCore.Qt.PenStyle.DashLine)))
        where = (f"undivided peak at {10 ** f.position:,.0f}"
                 + (f" (control: {10 ** self.peak.position:,.0f})" if self.peak else ", fitted"))
        self.said.setText(f"{self.shown}: {proliferation.summary(f)} · {where}")


# --- the cell cycle -------------------------------------------------------------------------

class CellCycleDialog(_ModelDialog):
    KIND, NAME, COMPARE = CELL_CYCLE, "cell_cycle", "% S"
    DEFINITIONS = (
        "G1 and G2 are Gaussians on the linear DNA signal. Watson pragmatic: each fitted on "
        "the side S does not reach; S is what is left between them, whatever its shape — "
        "noisier, and biased where S piles up against a peak. Dean-Jett-Fox: S is a "
        "broadened second-degree polynomial fitted with the peaks — the more precise where "
        "S is smooth, wrong where it is not. % G1, % S and % G2 are of the cells in the "
        "three; sub-G1 and above G2 are of all events, outside the dotted lines. Gate "
        "single cells first: two G1 cells stuck together count as G2. SE: standard error "
        "from counting.")

    def _controls(self) -> None:
        s = self.source
        self.model = QtWidgets.QComboBox()
        for key, text in cellcycle.MODELS.items():
            self.model.addItem(text, key)
        self.ratio_lo = _spin(cellcycle.RATIO[0], 1.5, 2.5, 2, 0.05)
        self.ratio_hi = _spin(cellcycle.RATIO[1], 1.5, 2.5, 2, 0.05)
        ratio = QtWidgets.QHBoxLayout()
        ratio.addWidget(self.ratio_lo)
        ratio.addWidget(QtWidgets.QLabel("to"))
        ratio.addWidget(self.ratio_hi)
        for box in (self.ratio_lo, self.ratio_hi):
            box.setToolTip("The mean of G2 is held between these multiples of the mean of "
                           "G1. The same number twice fixes it.")
        self.same_cv = QtWidgets.QCheckBox("G2 has the CV of G1")
        self.same_cv.setToolTip("Unticked, each peak has its own CV. Ticked, one CV for "
                                "both: steadier where G2 is small.")
        self.g1 = _spin(0.0, 0.0, 1e9, 0, 1000.0,
                        "About where the G1 peak is, in the units of the file, when the one "
                        "found is not it")
        self.g1.setSpecialValueText("Found")
        self.form.addRow("Model", self.model)
        self.form.addRow("G2 / G1", ratio)
        self.form.addRow("G1 near", self.g1)
        self.form.addRow(self.same_cv)
        self.model.currentIndexChanged.connect(self.changed)
        self.same_cv.toggled.connect(self.changed)
        for box in (self.ratio_lo, self.ratio_hi, self.g1):
            box.valueChanged.connect(self.changed)
        # DNA is read on a linear parameter: the one shown if it is one,
        # otherwise the first that is
        linear = [d for _, d in s.parameters if s.session.scale_of(d) is None
                  and d.strip().lower() != "time"]
        start = s.dim if s.dim in linear or not linear else linear[0]
        self.parameter.setCurrentIndex(max(self.parameter.findData(start), 0))

    def columns(self) -> list[str]:
        return list(cellcycle.COLUMNS)

    def set_model(self, kind: str) -> None:
        self.model.setCurrentIndex(max(self.model.findData(kind), 0))
        self.compute()

    def _fit(self, entry):
        f = cellcycle.fit(self.source.values(entry, self.parameter.currentData()),
                          self.model.currentData(),
                          (self.ratio_lo.value(), self.ratio_hi.value()),
                          self.same_cv.isChecked(), self.g1.value() or None)
        return f, cellcycle.row(f)

    def content(self):
        f = self.fits.get(self.shown)
        if f is None:
            return None
        return cellcycle.content(f, self.source.label(self.parameter.currentData()), self.shown)

    def _draw(self) -> None:
        content = self.content()
        if content is None:
            self._nothing()
            return
        self._draw_fit(content, linear_ticks=True)
        self.said.setText(f"{self.shown}: {cellcycle.summary(self.fits[self.shown])}")


# --- kinetics -------------------------------------------------------------------------------

class KineticsDialog(_ModelDialog):
    KIND, NAME, COMPARE = KINETICS, "kinetics", "Peak − baseline"
    DEFINITIONS = (
        "The acquisition is cut into slices of time; each point is the statistic of the "
        "population's events in one slice, the line is the points smoothed. Baseline: mean "
        "of the points in the grey band. Peak: the highest of the line in the response, "
        "between the dashed marks; Time to peak is counted from the first of them. Slope: "
        "of a straight line through the points from there to the peak. AUC: area between "
        "the line and the baseline over the response. Drag the bands or type their ends.")

    def _controls(self) -> None:
        s = self.source
        self._span: tuple | None = None
        self._unit = ""
        self.curves: list = []
        self.denominator = QtWidgets.QComboBox()
        self.denominator.addItem("—", "")
        for text, dim in s.parameters:
            self.denominator.addItem(text, dim)
        self.denominator.setToolTip("Follow the ratio of the parameter to this one, event "
                                    "by event: violet over blue of Indo-1, say")
        self.statistic = QtWidgets.QComboBox()
        for name in kinetics.STATISTICS:
            self.statistic.addItem(name, name)
        self.threshold = _spin(0.0, -1e9, 1e9, 2, 100.0,
                               "Events above this, in the units of the file (or of the "
                               "ratio), are the responding ones")
        self.slice = _spin(0.0, 0.0, 1e6, 2, 1.0,
                           "How long a slice of time is, in the units of the axis")
        self.slice.setSpecialValueText(f"{kinetics.SLICES} slices")
        self.smoothing = QtWidgets.QComboBox()
        for key, text in kinetics.SMOOTHINGS.items():
            self.smoothing.addItem(text, key)
        self.smoothing.setCurrentIndex(self.smoothing.findData(kinetics.MOVING))
        self.window = _spin(5.0, 1.0, 99.0, 0, 1.0,
                            "In slices: how many a moving average spans, or the SD of the "
                            "Gaussian")
        smooth = QtWidgets.QHBoxLayout()
        self.window.setFixedWidth(64)
        smooth.addWidget(self.smoothing, 1)
        smooth.addWidget(self.window)
        self.marks = [_spin(0.0, -1e9, 1e9, 2, 1.0) for _ in range(4)]
        rows = []
        for a, b in ((0, 1), (2, 3)):
            row = QtWidgets.QHBoxLayout()
            row.addWidget(self.marks[a])
            row.addWidget(QtWidgets.QLabel("to"))
            row.addWidget(self.marks[b])
            rows.append(row)
        self.form.addRow("Divided by", self.denominator)
        self.form.addRow("Statistic", self.statistic)
        self.form.addRow("Threshold", self.threshold)
        self.form.addRow("Slice", self.slice)
        self.form.addRow("Smoothing", smooth)
        self.form.addRow("Baseline", rows[0])
        self.form.addRow("Response", rows[1])
        self.threshold.setEnabled(False)
        self.statistic.currentIndexChanged.connect(self._statistic_changed)
        for combo in (self.denominator, self.smoothing):
            combo.currentIndexChanged.connect(self.changed)
        for box in (self.threshold, self.slice, self.window):
            box.valueChanged.connect(self.changed)
        for box in self.marks:
            box.valueChanged.connect(self._marks_typed)
        # the parameter shown, unless it is the time itself
        dims = [d for _, d in s.parameters if d.strip().lower() != "time"]
        start = s.dim if s.dim in dims or not dims else dims[0]
        self.parameter.setCurrentIndex(max(self.parameter.findData(start), 0))

    def _statistic_changed(self, *_) -> None:
        self.threshold.setEnabled(self.statistic.currentData() == kinetics.ABOVE)
        self.changed()

    def columns(self) -> list[str]:
        return list(kinetics.COLUMNS)

    def baseline(self) -> tuple[float, float]:
        return self.marks[0].value(), self.marks[1].value()

    def response(self) -> tuple[float, float]:
        return self.marks[2].value(), self.marks[3].value()

    def set_marks(self, baseline: tuple, response: tuple) -> None:
        self._filling = True
        for box, value in zip(self.marks, (*baseline, *response)):
            box.setValue(value)
        self._filling = False
        self._measure()

    def _marks_typed(self, *_) -> None:
        if not self._filling:
            self._measure()

    def ylabel(self) -> str:
        s = self.source
        what = s.label(self.parameter.currentData())
        if self.denominator.currentData():
            what += " / " + s.label(self.denominator.currentData())
        if self.statistic.currentData() == kinetics.ABOVE:
            return f"% of events with {what} above {self.threshold.value():g}"
        return f"{self.statistic.currentData()} of {what}"

    def _prepare(self) -> None:
        """Every sample's time, so that all are cut into the same slices."""
        self._times, steps, problems = {}, set(), {}
        for e in self.source.entries:
            try:
                self._times[e.name], step = self.source.time(e)
                steps.add(step)
            except ValueError as exc:
                problems[e.name] = str(exc)
        if not self._times:
            raise ValueError(next(iter(problems.values()), "no sample selected"))
        # in seconds only where every file says how long its unit is
        self._step = steps.pop() if len(steps) == 1 else 0.0
        self._unit = "s" if self._step else "as the file stores it"
        scale = self._step or 1.0
        have = [t for t in self._times.values() if t.size]
        if not have:
            raise ValueError("the population has no events")
        span = (min(float(t.min()) for t in have) * scale,
                max(float(t.max()) for t in have) * scale)
        if span != self._span:
            self._span = span
            self._filling = True
            for box, value in zip(self.marks, (v for pair in kinetics.default_marks(*span)
                                               for v in pair)):
                box.setSingleStep(max((span[1] - span[0]) / 100, 0.01))
                box.setValue(value)
            self._filling = False

    def _fit(self, entry):
        if entry.name not in self._times:
            self.source.time(entry)      # says why
        s = self.source
        bottom = self.denominator.currentData()
        c = kinetics.curve(
            self._times[entry.name], s.values(entry, self.parameter.currentData()),
            entry.name, s.group(entry), self._step, self.slice.value() or None,
            self.statistic.currentData(), self.threshold.value(),
            self.smoothing.currentData(), self.window.value(),
            s.values(entry, bottom) if bottom else None, self._span)
        return c, list(kinetics.measure(c, self.baseline(), self.response()).values())

    def _measure(self) -> None:
        """The marks moved: the same series, measured again."""
        if self.result is None:
            return
        first = self.result.columns.index("Path") + 1
        for row in self.result.rows:
            c = self.fits.get(row[0])
            if c is not None:
                row[first:first + len(kinetics.COLUMNS)] = list(
                    kinetics.measure(c, self.baseline(), self.response()).values())
        self._show_table(self.result)

    def content(self):
        curves = [self.fits[e.name] for e in self.source.entries if e.name in self.fits]
        if not curves:
            return None
        return kinetics.KineticsContent(curves, f"Time ({self._unit})", self.ylabel(),
                                        self.baseline(), self.response())

    def _draw(self) -> None:
        content = self.content()
        if content is None:
            self._nothing()
            return
        legend = self._legend()
        names = content.groups()
        colour = dict(zip(names, palette.categorical(len(names))))
        legend.setVisible(1 < len(names) <= 14)
        listed = set()
        for c in content.series:
            bold = c.label == self.shown or len(content.series) == 1
            dots = QtGui.QColor(colour[c.group])
            dots.setAlpha(90)
            self.plot.addItem(pg.ScatterPlotItem(c.t, c.y, size=4, pen=None,
                                                 brush=pg.mkBrush(dots)))
            line = pg.PlotCurveItem(c.t, c.smooth, connect="finite",
                                    pen=pg.mkPen(colour[c.group], width=2.4 if bold else 1.3))
            self.plot.addItem(line)
            if c.group not in listed:
                listed.add(c.group)
                legend.addItem(line, c.group)
        grey = QtGui.QColor(theme.muted())
        grey.setAlpha(45)
        clear = QtGui.QColor(theme.accent())
        clear.setAlpha(18)
        self.regions = []
        for k, (ends, brush) in enumerate(((self.baseline(), grey), (self.response(), clear))):
            region = pg.LinearRegionItem(ends, brush=pg.mkBrush(brush), pen=pg.mkPen(
                theme.muted(), width=1, style=QtCore.Qt.PenStyle.DashLine))
            region.setZValue(-10)
            region.sigRegionChangeFinished.connect(lambda r, k=k: self._dragged(k, r))
            self.plot.addItem(region)
            self.regions.append(region)
        self.plot.getAxis("bottom").setTicks(None)
        self.plot.setLabel("bottom", content.xlabel)
        self.plot.setLabel("left", content.ylabel)
        self.plot.enableAutoRange()
        shown = self.fits.get(self.shown)
        if shown is None:
            self._nothing()
            return
        m = kinetics.measure(shown, self.baseline(), self.response())
        self.said.setText(
            f"{self.shown}: baseline {m['Baseline']:.4g} · peak {m['Peak']:.4g} "
            f"({m['Peak / baseline']:.3g} × baseline) at {m['Time of peak']:.4g}, "
            f"{m['Time to peak']:.4g} after the response starts · slope {m['Slope']:.3g} · "
            f"AUC {m['AUC']:.4g}" + ("" if self._step else
                                     " · the files do not say how long a unit of Time is"))

    def _dragged(self, k: int, region) -> None:
        if self._filling:
            return
        lo, hi = region.getRegion()
        ends = [self.baseline(), self.response()]
        ends[k] = (lo, hi)
        self.set_marks(*ends)


DIALOGS = {PROLIFERATION: ProliferationDialog, CELL_CYCLE: CellCycleDialog,
           KINETICS: KineticsDialog}


def dialog_for(tab, kind: str):
    """The dialog of one model for what the Gates tab shows, fitted and
    not yet shown; None before there is a population to model."""
    if tab.sample is None or tab.node is None:
        tab.sigStatus.emit("Choose a sample first")
        return None
    dialog = DIALOGS[kind](Source(tab), tab)
    dialog.compute()
    return dialog
