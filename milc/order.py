"""The order of the values of a column: which group comes first on a chart,
in a legend, in a table.

Left alone, values come in the order a person counts in — day 3, 7, 14, not
14, 3, 7; `Case · 3` before `Case · 14`. The order of a column can also be
given (Ctrl before Case, or the other way round), and what is given wins.

The given order is kept here, for the whole program, rather than passed to
every function that sorts: there is one experiment open at a time, and its
order is the same on every chart. `Session` sets it when a folder is opened
and when it is changed. A value has one place whatever column it is in; if
two columns give the same value two places, the first given counts.
"""

from __future__ import annotations

import re

#: what joins the values of two columns into one label: `Case · 3`
JOIN = " · "

_rank: dict[str, int] = {}
_number = re.compile(r"(\d+(?:\.\d+)?)")


def set_order(orders: dict[str, list[str]] | None) -> None:
    """`orders[column]` is that column's values, first to last. Values not
    named come after those that are, in counting order."""
    _rank.clear()
    for values in (orders or {}).values():
        for k, value in enumerate(values):
            _rank.setdefault(str(value), k)


def natural(text: str) -> tuple:
    """A key under which `Day 3` comes before `Day 14`, and case does not
    matter."""
    parts = _number.split(str(text))
    return tuple((0, float(p), "") if k % 2 else (1, 0.0, p.casefold())
                 for k, p in enumerate(parts) if p != "")


def key(name) -> tuple:
    """What to sort by. A label made of several columns' values is ordered
    by the first, then by the second."""
    name = str(name)
    if name in _rank:
        return ((0, _rank[name], ()),)
    return tuple((0, _rank[p], ()) if p in _rank else (1, 0, natural(p))
                 for p in name.split(JOIN))


def ordered(names) -> list:
    """Distinct `names`, in order."""
    return sorted(dict.fromkeys(names), key=key)
