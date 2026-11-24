---
title: Keyboard shortcuts
---

MIL-C is mostly worked with the mouse: gates are drawn, samples clicked, values dragged. The
keyboard has a few shortcuts, listed here with the mouse gestures that are not obvious from the
buttons. Shortcuts follow each system's own conventions, so `Ctrl` on Windows and Linux is `⌘`
on macOS.

## The main window

| Action | Windows, Linux | macOS | Menu |
|---|---|---|---|
| Open a folder of `.fcs` files | `Ctrl+O` | `⌘O` | `File ▸ Open folder…` |
| Save the concatenated events | `Ctrl+S` | `⌘S` | `File ▸ Save concatenated file…` |
| The manual, on the page for where you are | `F1` | `⌘?` | `Help ▸ Manual` |

`F1` also works inside a dialog: pressed in the Quality dialog it opens [[quality-control]], in
the Figure dialog [[figures]], in the Compensation dialog [[compensation]]. Where a dialog has no
page of its own, the manual opens on the page of the tab behind it.

The other commands of the menus — `File ▸ Open workspace…`, `Open concatenated file…`, `Export
tables…`, `Export gates…`, `View ▸ Appearance`, `View ▸ Plate…`, `Help ▸ Export manual as
PDF…` — have no shortcut.

## The manual window

| Action | Windows, Linux | macOS |
|---|---|---|
| Search the manual | `Ctrl+F` | `⌘F` |
| Open the first search result | `Return` in the search box | `Return` |
| Back to the previous page | `Alt+←` | `⌘[` |
| Forward again | `Alt+→` | `⌘]` |

See [[welcome]] for how the search matches words.

## The Samples table

| Action | How |
|---|---|
| Select several rows | `Shift`-click for a run, `Ctrl`-click (`⌘`-click) for one more |
| Edit a cell of your own columns | double-click it, or select it and start typing; `F2` on Windows and Linux |
| Put one value into many samples | select their rows, then type into one of them |
| Cancel an edit | `Esc` |
| Rename a column of yours | double-click its header |
| Rename or remove a column | right-click its header |

See [[samples-tab]].

## Drawing gates

In the plot of the [[gates-tab]], once a gate tool is chosen:

| Tool | Gesture |
|---|---|
| **Rectangle** | drag across the events; on a histogram, a range |
| **Polygon** | click each corner; double-click to close |
| **Lasso** | hold the mouse button down and draw around the events; let go to close |
| **Ellipse** | drag the box it fits in, then turn and stretch it by its handles |
| **Quadrant** | click where the two lines should cross |
| any tool | `Esc` stops drawing |

The samples list of the Gates tab takes `Shift`-click and `Ctrl`-click (`⌘`-click) to select
several samples for an overlay. More in [[drawing-gates]].

## Maps and dialogs

- On the map of the [[clusters-tab]], `Shift`+drag outlines any set of events, freehand, and the
  tab reports the same numbers for it as for a cluster.
- In the Figure dialog, double-click a category's colour to change it; type over its label to
  rename it.
- In the Order of the groups dialog, drag a value to move it ([[sample-order]]).
- Dialogs close with `Esc` and accept with `Return` where they have an **OK** button, as everywhere
  on the system.
