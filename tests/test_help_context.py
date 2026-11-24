"""
F1 opens the page for where the user is.

A widget names its page once and everything inside it inherits the answer;
the shell resolves the focused widget, and a dialog resolves itself. What
is tested is the walk, the fallbacks, and that every page named by a
widget or a dialog is a page the manual has.
"""

import os
import re
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6 import QtCore, QtGui, QtWidgets  # noqa: E402

from milc.manual import manual  # noqa: E402
from milc.ui.help_window import HELP_PROPERTY, describe, help_page_for  # noqa: E402

UI = Path(__file__).resolve().parent.parent / "milc" / "ui"


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_the_nearest_named_ancestor_answers(qapp):
    outer = QtWidgets.QWidget()
    describe(outer, "gates-tab")
    middle = QtWidgets.QWidget(outer)
    inner = QtWidgets.QLineEdit(middle)
    assert help_page_for(inner) == "gates-tab"
    describe(middle, "compensation")
    assert help_page_for(inner) == "compensation"
    assert help_page_for(QtWidgets.QWidget()) is None
    assert help_page_for(None) is None
    assert inner.property(HELP_PROPERTY) is None


def test_every_page_a_dialog_names_is_in_the_manual():
    """Read from the source: some dialogs need a sample to be built at all."""
    named = set()
    for path in UI.glob("*.py"):
        text = path.read_text(encoding="utf-8")
        named |= set(re.findall(r'describe\(self, "([^"]+)"\)', text))
        named |= set(re.findall(r'PAGE = "([^"]+)"', text))
    assert len(named) >= 15
    missing = sorted(page for page in named if page not in manual().pages)
    assert missing == [], f"dialogs name pages the manual lacks: {missing}"


def test_the_shell_names_its_tabs_and_f1_opens_their_page(qapp):
    from milc.ui.shell import MainShell

    shell = MainShell()
    try:
        pages = manual().pages
        for tab, page in ((shell.samples, "samples-tab"), (shell.gates, "gates-tab"),
                          (shell.statistics, "statistics-tab"),
                          (shell.explore, "explore-tab"), (shell.clusters, "clusters-tab")):
            assert help_page_for(tab) == page and page in pages
        shell.tabs.setCurrentWidget(shell.gates)
        assert shell.context_page() == "gates-tab"
        event = QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_F1,
                                QtCore.Qt.KeyboardModifier.NoModifier)
        if event.matches(QtGui.QKeySequence.StandardKey.HelpContents):
            assert shell.eventFilter(shell.statistics, event) is True
            assert shell._help_window.current_page == "statistics-tab"
        shell.show_manual("compensation")
        assert shell._help_window.current_page == "compensation"
    finally:
        shell.close()
