---
title: Balanced concatenation
---

UMAP, tSNE and the clustering methods run on one set of events drawn from every sample: the
concatenated events. How many each sample gives decides what the map is a map of. A treated
sample often carries several times the cells of its control, and an embedding built on
everything is an embedding of the treated group. In MIL-C the plan — what each sample has, what
it will give, and which one set the limit — is on the table before anything is drawn, so the cost
of balancing is read before it is paid.

## The settings

In the **Downsampling** section of the [[samples-tab]]:

| Control | Meaning |
|---|---|
| **Balance** | one of three ways to decide each sample's number (below) |
| **N** | the number asked for; `as many as the data allows` (0) lets the data set it |
| **Minimum** | a sample with fewer events of the population than this is left out, instead of pulling every other sample down to its size |
| **Seed** | the seed of the random draw; the same seed gives the same events |

The events are those of the **Population** chosen above, in the samples ticked, grouped by the
columns of **Group by**.

## The three ways to balance

**Same number from every sample.** Every sample gives the same number: that of the smallest
sample, or **N** if it is smaller. The smallest sample's Note reads `sets the limit`. If **N** is
larger than the smallest sample has, the note `asked for N per sample; sample has only k` says
so and the smallest still sets the number.

**Same total from every group.** Every group gives the same total: that of the smallest group,
or **N** if it is smaller (here **N** is per group). Inside a group the total is spread as evenly
as the samples allow: a small sample gives everything it has, and the larger ones share what is
left, the odd events going to the largest samples. When samples of one group end up giving
different numbers the note says `inside a group the samples do not give the same number: the
larger ones make up for the smaller`.

**At most N from each sample.** Each sample gives **N**, or all it has if that is fewer. This
mode needs a number in **N**. It is not balanced when some samples are short, and says so:
`k sample(s) have fewer than N and give all they have: the samples are not balanced`.

## An example

Three Control and three Treated samples, with these events in the population:

| Sample | Group | Has | Same number | Same total per group | At most 20,000 |
|---|---|---|---|---|---|
| sample_01 | Control | 12,000 | 12,000 | 12,000 | 12,000 |
| sample_02 | Control | 15,000 | 12,000 | 15,000 | 15,000 |
| sample_03 | Control | 20,000 | 12,000 | 20,000 | 20,000 |
| sample_04 | Treated | 40,000 | 12,000 | 15,666 | 20,000 |
| sample_05 | Treated | 55,000 | 12,000 | 15,667 | 20,000 |
| sample_06 | Treated | 60,000 | 12,000 | 15,667 | 20,000 |
| **total** | | 202,000 | 72,000 | 94,000 | 107,000 |

- Same number: limited by sample_01; 36,000 per group; 36% of what the samples have.
- Same total: limited by the group Control (47,000); every Treated sample gives about a third
  of that.
- At most 20,000: the groups give 47,000 and 60,000, and the note warns that two samples give
  all they have.

Which to use: **Same number** when every sample should weigh the same — each animal is one
observation. **Same total** when the groups should weigh the same and samples are uneven within
them. **At most N** to cap very large files without discarding small ones, accepting that the
result is not balanced. **Minimum** keeps one nearly empty sample from cutting every other down:
with **Minimum** 5,000, a sample of 800 events is excluded and its Note says why.

## What the draw is

Each sample's events are drawn at random, without replacement, from the population, and kept in
the order they were acquired. The draw depends on the seed and on the sample's own name — not on
which other samples take part nor in what order — so unticking one sample does not change the
events drawn from the others.

## What Concatenate builds

The concatenated set holds, for every event drawn:

- every channel of the files, fluorescence compensated (with the matrix each sample is given),
  nothing transformed: the scale is chosen where the values are used;
- which sample it came from and which event of that file it was;
- the sample's annotation (every column of the table).

With it is kept how it was built: the population, the mode, what each sample had and gave, which
were excluded and why, which sample or group set the limit, and the seed. Where quality control
was applied ([[quality-control]]) the record says what was removed from each sample; where batch
normalisation is on ([[normalisation]]) the batches are aligned before the set is finished, and
the record says by how much. All of this appears in the report of [[the-plate]] and is saved
with the set ([[projects-and-saving]]).

After **Concatenate** the window moves to [[explore-tab]]. A new concatenation replaces the
previous set and the results computed on it: save the set first if you want to keep it.

## When it goes out of date

If the values under the samples change after the set was built — another compensation matrix,
unmixing, normalisation applied to the samples themselves — a line above the tabs says that
Explore and Clusters still show the old values until you concatenate again. Changing the plan
itself (ticks, **N**, the population) does not touch a set already built.

## From the command line

`milc concat` does the same without the window, prints the plan, and writes `NAME.npz` and
`NAME.fcs`; `--dry-run` prints the plan only:

```
milc concat DATA --wsp analysis.wsp --population Cells/Live/CD45 \
    --mode equal_per_group --group-by Group --min-events 1000 --seed 0 --out concat
```

The modes are `equal_per_sample`, `equal_per_group` and `ceiling`. See [[command-line]].
