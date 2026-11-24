---
title: The order of the groups
---

Which group comes first is a small thing that shows everywhere: the bars of a frequency chart,
the entries of a legend, the panels of a split map, the rows of a table, the columns of a Prism
file. MIL-C keeps one order for each column of the Samples table and uses it on every one of
those, so that Control stands left of Treated on the map, the chart and the exported table
alike.

## Counting order

Left alone, values come in the order a person counts in, not the order a computer sorts text:

| Values | Counting order | Plain text order |
|---|---|---|
| day 3, 7, 14 | 3, 7, 14 | 14, 3, 7 |
| `sample_2`, `sample_10` | `sample_2`, `sample_10` | `sample_10`, `sample_2` |
| Treated, control | control, Treated: case is ignored | Treated, control: capitals first |

Numbers inside a value are compared as numbers, and capital and small letters as the same
letter. The files of a folder are listed in the same order.

## Giving an order

**Order…** opens the dialog **Order of the groups**. It is in the Samples tab, above the table,
and in the Clusters tab, where it is most often wanted.

- **Column** chooses the column whose values are shown. The dialog opens on the first column, or
  on the one the button was pressed for.
- Drag a value up or down the list, or select it and use **Move up** and **Move down**.
- **Reverse** turns the list upside down.
- **Counting order** goes back to the order nobody chose: 3, 7, 14; A, B, C.
- The line under the list says which it is: `In the order given here.`, `In counting order:
  nobody chose another.`, or `No values in this column yet.`

Each column keeps its own order. You can arrange one column, switch to another and arrange that,
and everything arranged is applied together on **OK**. **Cancel** changes nothing.

## Where it applies

Once applied, the order is used by every part of the program that shows groups:

- the frequency charts, violins and boxes, and the legends of every plot;
- the panels of a map split by a column, in [[explore-tab]] and [[clusters-tab]];
- the groups of the tests between groups ([[statistics-tab]], [[cluster-statistics]]);
- the tables of a Prism file ([[tables-and-prism]]);
- the curves of the batch-normalisation dialog ([[normalisation]]) and the choices of
  **Reference** and **Bring to**;
- figures and plates made from any of these ([[figures]]).

Nothing is computed again when the order changes: the order is of what is shown, not of what was
measured. Every plot that shows groups is simply drawn again, and the status bar says `The order
of the groups changed`.

## Groups made of two columns

When **Group by** ticks two columns, a group is a combination, written with a dot: `Treated ·
7`. Such a label is ordered by its first column, then by its second, each in its own order. With
Group in the order Control, Treated and Day in counting order, the groups come as `Control · 3`,
`Control · 7`, `Control · 14`, `Treated · 3`… Which column is outermost follows the order of the
columns in the Samples table, whatever order they were ticked in.

## Rules worth knowing

- A value has one place, whatever column it is in. If two columns contain the same value and
  give it different places, the first column's order counts.
- Values not named in a column's order — a new group typed into the table later — come after the
  named ones, in counting order.
- Renaming a column keeps its order; removing a column removes its order.
- The order is kept in the folder's project file ([[projects-and-saving]]) and is there next time
  the folder is opened.
