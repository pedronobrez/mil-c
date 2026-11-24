---
title: Quality control of the acquisition
---

A clog, an air bubble, a change of pressure or a drifting laser spoils the events acquired while
it lasts. MIL-C looks for such stretches the way flowAI and PeacoQC do — in the order the events
were acquired, a fixed number of events at a time — and what it finds is a list of time
intervals. Applied, they become a population, **Clean**, at the top of every sample's gating
tree: every gate, statistic and concatenation then starts from the clean events. No R is needed.

## Opening it

**Quality…** is in the [[gates-tab]], among the gate tools; choose a sample first. The dialog
**Quality of the acquisition: sample** opens on that sample, already checked with the limits last
applied (or the defaults). Nothing is removed by looking.

## What is checked

The events are put in the order of the **Time** parameter and cut into bins of the same number
of events. A bin is suspect when:

| Check | What it finds |
|---|---|
| **Flow rate** | the bin's events per second are far from the acquisition's usual rate: a clog, a bubble. A rate within 1.5 times the usual one is never a defect |
| **Signal** | the running median, over five bins, of a population's bin medians lies more than the MAD limit away from that population's median over the whole file: a step, anything abrupt |
| **Drift** | the slow level of a population, followed over a long window, moves away from where it stays for longest: a moderate shift that lasts |
| **Shift** | a short stretch, a few bins long, whose level stands apart from the bins on both sides: a moderate shift that does not last |

Then a stretch of fewer than **Consecutive bins** good bins between two removed stretches is
removed too.

Every parameter but Time is judged, fluorescence compensated and each on the scale its axis is
drawn on. A parameter with several populations — a marker half the cells carry — is cut at the
valleys of its histogram and followed one population at a time, so that a bin with a slightly
different mix of positive and negative cells is not taken for a defect. A population too small
to give a median in most bins (fewer than 20 events in half of them) is not followed. A defect is
found only while it touches less than half of the acquisition: with more, the bad part is the
reference.

## The limits

| Control | Default | Meaning |
|---|---|---|
| **MAD limit** | 6 | how many MADs from the acquisition's own level the signal, or the flow rate, has to lie for a bin to be suspect. Larger removes less |
| **Consecutive bins** | 5 | good bins between two removed stretches are kept only when there are at least this many |
| **Events per bin** | auto | auto is 1/500 of the file and at least 150 events. Larger bins are steadier and see less of a short defect |
| **Drift limit** | 0.25 | the smallest lasting shift taken for a drift, as a share of the spread of the parameter's events. **off** (0) does not look for drifts |
| **Shift limit** | 0.25 | the smallest move of a short stretch taken for a shift, in the same unit; it must also stand out by the MAD limit. **off** does not look |

Change a limit and **Check again** checks this sample with the new limits — and every sample, if
they were checked with the old ones.

## Reading the plots

The line at the top summarises what would go: `To remove: 3.2% of 120,000 events, in 2
stretch(es) of time — flow rate 1,200 events; signal 2,640 (CD4, CD8).` Where more than a fifth
of the file would go it adds: `look at the plots before applying.`

- **Top plot**: the flow rate of each bin against time (in seconds when the file gives
  `$TIMESTEP`, otherwise in the file's own units). Dashed lines: how far it may go; red points:
  bins beyond.
- **Bottom plot**, for the population chosen in **Signal of**: the median of each bin, with its
  dashed limits and red suspect bins; a thick line, the slow level, red where it drifted beyond
  the dotted limits; red bars over short stretches that stood apart from their neighbours (dashed
  bars).
- **Shaded** in both: the stretches to remove.

**Signal of** starts on the population that went wrong for longest. The time axes of the two
plots move together.

## Every sample

**Check every sample** checks each file with the current limits, off the window's thread, and
fills a table: **Sample**, **Events**, **To remove** (%), the events found by **Flow rate**,
**Signal**, **Drift** and **Shift**, and a **Note** naming the parameters to blame — or why the
file could not be checked. A sample that could be checked has a box, ticked.

A file cannot be checked when it has no Time parameter, has too few events (it takes 20 bins:
3,000 at the smallest bin), has Time missing for some events, has the same Time for every event,
or has Time running backwards (several acquisitions in one file, or a clock that was reset).

## Applying and removing

- **Apply to this sample** removes the stretches found in the sample shown.
- **Apply to every sample** (after **Check every sample**) removes from each ticked sample its own
  stretches; an unticked sample keeps every event.
- **Remove** takes **Clean** out of every sample's tree; the gates go back to where they were.

Applying puts **Clean** at the top of every sample's tree — samples left whole get it too, keeping
all their events, so the trees stay alike — and every existing gate moves under it. The
population chosen in the Samples tab moves with it (`Cells/Live` becomes `Clean/Cells/Live`).
Automatic gates fitted on the data are fitted again on the clean events ([[auto-gates]]). The
line at the foot of the dialog says what the trees hold, and the limits last applied; these are
saved with the project and written into the report of [[the-plate]].

## How it compares

Run beside PeacoQC 1.22.0 and flowAI 1.42.0 on the same values (details in [[measured-facts]]):

- On known defects written into a clean file and on simulated ones, MIL-C removed a step and a
  clog completely; PeacoQC, which has no flow-rate check, missed the clogs.
- On 32 real files, MIL-C removed 5.25% of all events, PeacoQC 5.13% and flowAI 11.97%.
- The drift and shift limits were measured on copies of real files with shifts of known size: the
  real files lose the same with **Drift limit** anywhere from 0.20 to 0.33 and **Shift limit**
  from 0.20 to 0.30; 0.25 is the middle of both.

## From the command line

`milc qc DATA --out qc.csv --plots qc.pdf` checks every file and writes the table and a page of
plots per sample, removing nothing. `--qc` on `concat`, `stats`, `export` and `gates` removes the
stretches for that run only, without writing the project. See [[command-line]].
