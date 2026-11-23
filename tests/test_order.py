"""Which group comes first: counting order where nobody chose, the user's
order where they did, on charts, legends and tables alike."""

import numpy as np
import pytest

from milc import cluster_stats, compare, figures, order
from milc.session import Session
from test_ui import app, batches  # noqa: F401  (fixtures)

PATTERN = r"(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)"


@pytest.fixture(autouse=True)
def nothing_given():
    order.set_order(None)
    yield
    order.set_order(None)


def test_left_alone_values_come_as_a_person_counts():
    assert order.ordered(["14", "3", "7", "3"]) == ["3", "7", "14"]
    assert order.ordered(["Day 14", "day 3", "Day 7"]) == ["day 3", "Day 7", "Day 14"]
    assert order.ordered(["Ctrl", "Case"]) == ["Case", "Ctrl"]
    # two columns as one label: by the first, then by the second
    assert order.ordered(["Ctrl · 3", "Case · 14", "Case · 3", "Ctrl · 14"]) == [
        "Case · 3", "Case · 14", "Ctrl · 3", "Ctrl · 14"]
    assert order.ordered(["b", "", "a"]) == ["", "a", "b"]
    assert order.ordered([2, 10, 1]) == [1, 2, 10]


def test_an_order_that_was_given_wins_everywhere():
    order.set_order({"Group": ["Ctrl", "Case"], "Day": ["14", "3"]})
    assert order.ordered(["Case", "Ctrl"]) == ["Ctrl", "Case"]
    assert order.ordered(["3", "7", "14"]) == ["14", "3", "7"]        # 7 was not named: after
    assert order.ordered(["Case · 3", "Ctrl · 14", "Ctrl · 3", "Case · 14"]) == [
        "Ctrl · 14", "Ctrl · 3", "Case · 14", "Case · 3"]
    # the two functions the tables and the tests use are this one
    assert cluster_stats.ordered(["Case", "Ctrl"]) == compare.ordered(["Case", "Ctrl"]) == [
        "Ctrl", "Case"]
    # the legend of a map, and so its colours
    c = figures.categories_of(np.array(["Case", "Ctrl", "Case", "Ctrl"]))
    assert c.names == ["Ctrl", "Case"] and list(c.codes) == [1, 0, 1, 0]
    order.set_order(None)
    assert figures.categories_of(np.array(["Case", "Ctrl"])).names == ["Case", "Ctrl"]


def test_the_order_is_the_projects(batches):
    s = Session()
    s.open_folder(batches)
    s.annotate(PATTERN)
    assert s.values_of("Group") == ["Case", "Ctrl"] and s.values_of("Day") == ["3", "7"]
    s.set_order("Group", ["Ctrl", "Case"])
    s.set_order("Day", ["7", "3"])
    assert s.values_of("Group") == ["Ctrl", "Case"] and s.values_of("Day") == ["7", "3"]

    again = Session()
    again.open_folder(batches)
    assert again.orders == {"Group": ["Ctrl", "Case"], "Day": ["7", "3"]}
    assert again.values_of("Group") == ["Ctrl", "Case"]
    # a column renamed keeps its order; counting order is asked for with nothing
    again.rename_column("Day", "Time")
    assert again.orders["Time"] == ["7", "3"]
    again.set_order("Time", None)
    assert again.values_of("Time") == ["3", "7"] and "Time" not in again.orders
    # another folder starts with nothing given
    again.open_folder(str(batches))
    again.remove_column("Group")
    assert "Group" not in again.orders and order.ordered(["Ctrl", "Case"]) == ["Case", "Ctrl"]


def test_the_dialog_and_the_charts(app, batches):
    from milc.ui import worker
    from milc.ui.order_dialog import OrderDialog
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(4):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(batches)
    w.show()
    s = w.session
    w.samples.fill_from_names(PATTERN)
    w.samples.concatenate()
    settle()
    w.explore.clusters.setValue(3)
    w.explore.grid.setValue(4)
    w.explore.run(["FlowSOM"])
    settle()
    w.tabs.setCurrentWidget(w.clusters)
    settle()
    assert [n for n, _ in figures.category_names(w.clusters.freq.content())] == ["Case", "Ctrl"]

    d = OrderDialog(s, w, "Group")
    assert d.values() == ["Case", "Ctrl"] and "counting order" in d.state.text()
    d.list.setCurrentRow(1)
    d.move(-1)
    assert d.values() == ["Ctrl", "Case"] and "given here" in d.state.text()
    d.set_column("Day")
    assert d.values() == ["3", "7"]
    d.reverse_all()
    d.set_column("Group")                       # what was arranged for it is still there
    assert d.values() == ["Ctrl", "Case"]
    d.set_column("Day")
    d.counting_order()                          # and day goes back to nobody's order
    assert d.values() == ["3", "7"] and "counting order" in d.state.text()
    assert d.apply() and s.orders == {"Group": ["Ctrl", "Case"]}
    assert not OrderDialog(s, w, "Group").apply()         # nothing arranged: nothing changes

    w.order_changed()
    settle()
    assert [n for n, _ in figures.category_names(w.clusters.freq.content())] == ["Ctrl", "Case"]
    assert w.clusters.order_button.isVisible() and w.samples.order_button is not None
    w.close()
