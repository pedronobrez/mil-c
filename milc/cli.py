"""
The command line: everything the window will do, without the window.

    python -m milc                    the window
    python -m milc inspect  DATA
    python -m milc counts   DATA --wsp analysis.wsp
    python -m milc concat   DATA --wsp analysis.wsp --population PATH --out NAME
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

from . import __version__, cluster_stats, sampling
from .dataset import Dataset, concatenate
from .fcs import read_header
from .sample import Sample
from .transforms import Arcsinh, Biex

#: `14Case4` → day 14, group Case, animal 4; a trailing ` marks a re-run
DEFAULT_PATTERN = r"(?P<day>\d+)(?P<group>[A-Za-z]+)(?P<animal>\d+)(?P<rerun>`)?$"


def _files(folder: str) -> list[str]:
    return sorted(glob.glob(os.path.join(folder, "*.fcs")))


def cmd_inspect(a) -> int:
    for path in _files(a.data):
        kv = read_header(path)
        n = int(kv["$PAR"])
        stains = [kv.get(f"$P{i}S", "").strip() for i in range(1, n + 1)]
        print(f"{os.path.basename(path):30s} {int(kv['$TOT']):>9,d} events  "
              f"{n} parameters  {kv.get('$CYT', '')}  {kv.get('$DATE', '')}")
        if a.verbose:
            for i in range(1, n + 1):
                print(f"    {kv[f'$P{i}N']:18s} {stains[i - 1]}")
    return 0


def _load(a) -> tuple[list[Sample], list[str]]:
    """The samples that can be gated, through the same session the window
    uses; and the files of the folder that cannot."""
    from .session import Session
    session = Session()
    session.open_folder(a.data)
    session.set_workspace(a.wsp)
    session.set_extend_gates(not a.listed_only)
    session.annotate(a.pattern)
    samples, left_out = [], []
    for e in session.entries:
        if e.gated and not e.problem:
            s = session.load(e)
            s.meta = dict(e.meta)
            samples.append(s)
        else:
            left_out.append(os.path.basename(e.path))
    return samples, left_out


def cmd_counts(a) -> int:
    samples, unknown = _load(a)
    worst = 0.0
    for s in samples:
        print(f"{s.name}  ({len(s):,d} events)")
        for node in s.gates.walk():
            if node.gate is None:
                continue
            ours, ref = s.count(node.path), node.reference_count
            depth = node.path.count("/")
            line = f"  {'  ' * depth}{node.name:<{22 - 2 * depth}s} {ours:>8,d}"
            if ref is not None:
                rel = (ours - ref) / ref if ref else 0.0
                worst = max(worst, abs(rel))
                line += f"   FlowJo {ref:>8,d}   {ours - ref:+d} ({rel:+.2%})"
            print(line)
    print(f"\nlargest difference from FlowJo: {worst:.2%}")
    if unknown:
        print(f"{len(unknown)} file(s) in the folder are not gated: "
              f"{', '.join(unknown)}")
    return 0


def cmd_concat(a) -> int:
    samples, _ = _load(a)
    if a.population is None:
        paths = {n.path for s in samples for n in s.gates.leaves()}
        if len(paths) != 1:
            print("choose --population; the samples end in: "
                  + ", ".join(sorted(paths)), file=sys.stderr)
            return 2
        a.population = paths.pop()
    available = {s.name: s.count(a.population) for s in samples}
    groups = {s.name: "/".join(s.meta.get(k, "") for k in a.group_by.split(","))
              for s in samples}
    p = sampling.plan(available, groups, a.mode, a.n, a.min_events)

    print(f"population: {a.population}\nmode: {p.mode}"
          + (f", limited by {p.limited_by}" if p.limited_by else ""))
    for s in samples:
        why = f"   excluded: {p.excluded[s.name]}" if s.name in p.excluded else ""
        print(f"  {s.name:12s} {groups[s.name]:12s} {available[s.name]:>8,d} "
              f"→ {p.take[s.name]:>7,d}{why}")
    for g, t in sorted(p.group_totals().items()):
        print(f"  group {g:12s} {t:>8,d}")
    print(f"total {p.total:,d} events, {p.used_fraction():.1%} of what the "
          "included samples have")
    for note in p.notes:
        print(f"note: {note}")
    if a.dry_run:
        return 0
    ds = concatenate(samples, p, a.population, a.seed)
    ds.save(a.out + ".npz")
    ds.to_fcs(a.out + ".fcs")
    print(f"wrote {a.out}.npz and {a.out}.fcs")
    return 0


def cmd_analyze(a) -> int:
    from . import embed, flowsom  # numba and friends, only when asked for

    ds = (Dataset.load(a.dataset) if a.dataset.endswith(".npz")
          else Dataset.from_fcs(a.dataset, a.sample_column))
    if not a.channels:
        print("choose --channels among:", file=sys.stderr)
        for c, m in zip(ds.channels, ds.markers):
            print(f"  {c:18s} {m}", file=sys.stderr)
        return 2
    channels = [ds.resolve(c) for c in a.channels.split(",")]
    if a.transform == "biex":
        transform, scale = Biex(), 1023.0
    else:
        transform, scale = Arcsinh(a.cofactor), 1.0
    X = ds.matrix(channels, transform, scale)
    base = {"channels": channels, "transform": a.transform, "seed": a.seed,
            **({"cofactor": a.cofactor} if a.transform == "arcsinh" else {})}
    print(f"{len(ds):,d} events, {len(ds.samples)} samples; on "
          + ", ".join(ds.label(c) for c in channels) + f" ({a.transform})")

    if a.umap:
        Y = embed.umap(X, a.neighbors, a.min_dist, seed=a.seed)
        ds.derived["UMAP1"], ds.derived["UMAP2"] = Y[:, 0], Y[:, 1]
        ds.runs["UMAP"] = {**base, "n_neighbors": a.neighbors, "min_dist": a.min_dist}
        print(f"UMAP: keeps {embed.neighbour_preservation(X, Y):.1%} of each "
              "event's 15 nearest neighbours")
    if a.tsne:
        Y = embed.tsne(X, a.perplexity, seed=a.seed)
        ds.derived["tSNE1"], ds.derived["tSNE2"] = Y[:, 0], Y[:, 1]
        ds.runs["tSNE"] = {**base, "perplexity": a.perplexity}
        print(f"tSNE: keeps {embed.neighbour_preservation(X, Y):.1%} of each "
              "event's 15 nearest neighbours")
    if a.flowsom:
        fs = flowsom.flowsom(X, a.clusters, a.grid, a.grid, seed=a.seed)
        ds.derived["FlowSOM"] = fs.labels.astype(np.int32)
        ds.derived["FlowSOM_node"] = fs.node.astype(np.int32)
        ds.runs["FlowSOM"] = {**base, "n_clusters": a.clusters, "grid": a.grid}
        labels = fs.labels
        sizes = np.bincount(labels)[1:]
        print("FlowSOM: " + ", ".join(
            f"{i + 1}: {n:,d} ({n / len(ds):.1%})" for i, n in enumerate(sizes)))
        names = [ds.label(c) for c in channels]
        cluster_stats.expression(X, labels, names).to_csv(a.out + "_mfi.csv", "cluster")
        freq = cluster_stats.frequencies(labels, ds.sample_index, ds.samples)
        freq.to_csv(a.out + "_freq.csv", "cluster")
        groups = ["/".join(m.get(k, "") for k in a.group_by.split(","))
                  for m in ds.meta]
        if len(set(groups)) > 1:
            mean, sd, n = cluster_stats.by_group(freq, groups)
            mean.to_csv(a.out + "_freq_mean.csv", "cluster")
            sd.to_csv(a.out + "_freq_sd.csv", "cluster")

    ds.save(a.out + ".npz")
    keys = list(ds.derived)
    with open(a.out + "_events.csv", "w", encoding="utf-8") as f:
        f.write(",".join(["sample", "event"] + keys) + "\n")
        cols = [ds.derived[k] for k in keys]
        for i in range(len(ds)):
            f.write(",".join([ds.samples[ds.sample_index[i]], str(ds.event_index[i] + 1)]
                             + [f"{c[i]:.6g}" for c in cols]) + "\n")
    print(f"wrote {a.out}.npz, {a.out}_events.csv"
          + (" and the cluster tables" if a.flowsom else ""))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="milc")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="command", required=True)

    s = sub.add_parser("inspect", help="what a folder of FCS files holds")
    s.add_argument("data")
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(run=cmd_inspect)

    def common(s):
        s.add_argument("data")
        s.add_argument("--wsp", required=True)
        s.add_argument("--listed-only", action="store_true",
                       help="only the files the workspace lists; by default the "
                            "others are given the gates and the matrix those share")
        s.add_argument("--pattern", default=DEFAULT_PATTERN,
                       help="regular expression with named groups, read "
                            "from each file name")

    s = sub.add_parser("counts", help="every gate, against FlowJo's count")
    common(s)
    s.set_defaults(run=cmd_counts)

    s = sub.add_parser("concat", help="balanced downsampling and concatenation")
    common(s)
    s.add_argument("--population", help="path of the gate, e.g. Cells/Live/CD45")
    s.add_argument("--mode", choices=sampling.MODES, default=sampling.EQUAL_PER_SAMPLE)
    s.add_argument("--n", type=int)
    s.add_argument("--min-events", type=int, default=1)
    s.add_argument("--group-by", default="group",
                   help="annotation(s) that define a group, comma-separated")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--dry-run", action="store_true", help="the plan only")
    s.add_argument("--out", default="concat")
    s.set_defaults(run=cmd_concat)

    s = sub.add_parser("analyze", help="UMAP, tSNE and FlowSOM on a concatenated file")
    s.add_argument("dataset", help="the .npz from concat, or any concatenated .fcs")
    s.add_argument("--channels", help="detectors or stains, comma-separated")
    s.add_argument("--transform", choices=("biex", "arcsinh"), default="biex",
                   help="biex is what FlowJo feeds its plugins")
    s.add_argument("--cofactor", type=float, default=150.0)
    s.add_argument("--umap", action="store_true")
    s.add_argument("--neighbors", type=int, default=15)
    s.add_argument("--min-dist", type=float, default=0.5)
    s.add_argument("--tsne", action="store_true")
    s.add_argument("--perplexity", type=float, default=30.0)
    s.add_argument("--flowsom", action="store_true")
    s.add_argument("--clusters", type=int, default=5)
    s.add_argument("--grid", type=int, default=10)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--group-by", default="group")
    s.add_argument("--sample-column", default="SampleID",
                   help="for a .fcs: the column that says which sample")
    s.add_argument("--out", default="analysis")
    s.set_defaults(run=cmd_analyze)

    s = sub.add_parser("gui", help="the window (also what no command at all opens)")
    s.add_argument("folder", nargs="?")
    s.set_defaults(run=lambda a: __import__("milc.app", fromlist=["main"]).main(a.folder))

    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        argv = ["gui"]
    a = ap.parse_args(argv)
    return a.run(a)


if __name__ == "__main__":
    sys.exit(main())
