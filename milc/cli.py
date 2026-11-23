"""
The command line: everything the window will do, without the window.

    python -m milc                    the window
    python -m milc inspect  DATA
    python -m milc counts   DATA --wsp analysis.wsp
    python -m milc concat   DATA --wsp analysis.wsp --population PATH --out NAME
    python -m milc qc       DATA --out qc.csv --plots qc.pdf
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
from .transforms import Arcsinh, Biex, Logicle



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


def _load(a) -> tuple[list[Sample], list[str], object]:
    """The samples that can be gated, through the same session the window
    uses; the files of the folder that cannot; and the session."""
    from .session import Session
    session = Session()
    session.open_folder(a.data)
    session.set_workspace(a.wsp)
    session.set_extend_gates(not a.listed_only)
    if a.pattern:
        session.annotate(a.pattern)
    gated = _gated(session)
    reasons = _clean(session, a)
    samples, left_out = [], []
    for e in session.entries:
        if e in gated:
            s = session.load(e)
            s.meta = dict(e.meta)
            samples.append(s)
        else:
            left_out.append(os.path.basename(e.path))
    _say_removed(session.quality_record(samples), reasons)
    return samples, left_out, session


def cmd_counts(a) -> int:
    try:
        samples, unknown, _ = _load(a)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
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
    try:
        samples, _, session = _load(a)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    a.population = _under_clean(session, a.population)
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
    if a.normalize_by:
        column, _, value = (a.reference or "").partition("=")
        if not value:
            print("--normalize-by needs --reference COLUMN=VALUE: the samples that "
                  "should look alike in every batch", file=sys.stderr)
            return 2

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
    cleaned = session.quality_record([s for s in samples if p.take.get(s.name, 0) > 0])
    if cleaned:
        ds.provenance["quality"] = cleaned
    if a.normalize_by:
        # the session the samples came from: its gates are the ones in use,
        # the population of quality control among them
        s = session
        s.population, s.seed = a.population, a.seed
        s.normalization = {"batch": a.normalize_by, "reference": [column, value],
                           "clusters": a.normalize_clusters}
        if a.normalize_to:
            s.normalization["goal"] = a.normalize_to
        try:
            s._normalise(ds, {})
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            return 2
        done = ds.provenance["normalization"]
        print(f"normalised by {a.normalize_by} on {column}={value}; references apart by "
              f"{np.mean(done['distance_before']):.4f} before and "
              f"{np.mean(done['distance_after']):.4f} after (mean over markers, scale units)")
    ds.save(a.out + ".npz")
    ds.to_fcs(a.out + ".fcs")
    print(f"wrote {a.out}.npz and {a.out}.fcs")
    if cleaned:
        _quality_csv(a.out + "_qc.csv", [
            _quality_row(name, v["events"], v["checked"] is not False, v["removed"],
                         v["intervals"], v["timestep"])
            for name, v in cleaned["samples"].items()])
        print(f"wrote {a.out}_qc.csv: what quality control removed from each sample")
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
    elif a.transform == "logicle":
        transform, scale = Logicle(), 1023.0
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

    if a.phenograph:
        from . import graphcluster
        found = graphcluster.phenograph(X, a.k, a.resolution, a.seed, method=a.method,
                                        min_cluster_size=a.min_cluster_size)
        ds.derived["PhenoGraph"] = found.labels.astype(np.int32)
        ds.runs["PhenoGraph"] = {**base, **found.settings()}
        print(f"PhenoGraph: {found.n_clusters} clusters, modularity "
              f"{found.modularity:.3f}; the largest holds "
              f"{np.mean(found.labels == 1):.1%} of the events"
              f" ({graphcluster.describe(ds.runs['PhenoGraph'])})")
        _cluster_tables(ds, X, found.labels, channels, a.group_by,
                        a.out + "_phenograph")
    if a.pca:
        found = embed.pca(X, 2)
        ds.derived["PCA1"], ds.derived["PCA2"] = found.scores[:, 0], found.scores[:, 1]
        ds.runs["PCA"] = {"channels": channels, "transform": a.transform,
                          "explained": [float(v) for v in found.explained],
                          "loadings": [[float(v) for v in row] for row in found.loadings]}
        print("PCA: " + ", ".join(f"PCA{i + 1} holds {v:.1%} of the variance"
                                  for i, v in enumerate(found.explained)))
        cluster_stats.Table(
            [f"PCA{i + 1}" for i in range(len(found.loadings))],
            [ds.label(c) for c in channels], found.loadings,
        ).to_csv(a.out + "_pca_loadings.csv", "component")

    ds.save(a.out + ".npz")
    keys = list(ds.derived)
    with open(a.out + "_events.csv", "w", encoding=cluster_stats.CSV_ENCODING) as f:
        f.write(",".join(["sample", "event"] + keys) + "\n")
        cols = [ds.derived[k] for k in keys]
        for i in range(len(ds)):
            f.write(",".join([ds.samples[ds.sample_index[i]], str(ds.event_index[i] + 1)]
                             + [f"{c[i]:.6g}" for c in cols]) + "\n")
    print(f"wrote {a.out}.npz, {a.out}_events.csv"
          + (" and the cluster tables" if a.flowsom or a.phenograph else ""))
    return 0


def _cluster_tables(ds, X, labels, channels, group_by: str, stem: str) -> None:
    """What FlowSOM's clusters are written as above, for another
    clustering: medians, frequencies per sample, and their mean and SD
    per group when there are groups."""
    names = [ds.label(c) for c in channels]
    cluster_stats.expression(X, labels, names).to_csv(stem + "_mfi.csv", "cluster")
    freq = cluster_stats.frequencies(labels, ds.sample_index, ds.samples)
    freq.to_csv(stem + "_freq.csv", "cluster")
    groups = ["/".join(m.get(k, "") for k in group_by.split(",")) for m in ds.meta]
    if len(set(groups)) > 1:
        mean, sd, _ = cluster_stats.by_group(freq, groups)
        mean.to_csv(stem + "_freq_mean.csv", "cluster")
        sd.to_csv(stem + "_freq_sd.csv", "cluster")


def cmd_selftest(a) -> int:
    """The whole route on invented data, with no screen: what a build is
    asked before it is shipped. Says what it did and how long it took.

    With `--log FILE` everything goes to the file as well and a failure is
    written there and returned as an exit code, never raised: a windowed
    program on Windows has no console, and an error raised in one is a
    dialog that a build machine cannot dismiss."""
    if not a.log:
        return _selftest(a, print)
    import traceback
    with open(a.log, "w", encoding="utf-8") as log:
        def say(*parts, **_):
            text = " ".join(str(p) for p in parts)
            log.write(text + "\n")
            log.flush()
            print(text, flush=True)
        try:
            return _selftest(a, say)
        except BaseException:
            log.write(traceback.format_exc())
            return 1


def _selftest(a, print) -> int:
    import tempfile
    import time

    from . import figures
    from .fcs import write_fcs
    from .session import Session

    started = time.time()

    def step(text):
        print(f"  {time.time() - started:5.1f}s  {text}", flush=True)

    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    with tempfile.TemporaryDirectory() as folder:
        for name, n in (("3Ctrl1", 300), ("3Ctrl2", 450), ("3Case1", 900), ("3Case2", 1200)):
            pop = rng.integers(0, 3, n)
            fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
            ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
            write_fcs(os.path.join(folder, name + ".fcs"), ev,
                      ["FSC-A", "FITC-A", "PE-A", "APC-A"], markers=["", "CD8", "CD4", "Foxp3"],
                      keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
        s = Session()
        s.open_folder(folder)
        ds = s.build()
        assert len(ds) == 1200 and len(ds.samples) == 4, ds.counts()
        step(f"read 4 files, concatenated {len(ds):,d} events")
        channels = s.default_channels()
        s.run_umap(channels)
        step(f"UMAP, keeps {ds.runs['UMAP']['neighbours_kept']:.0%} of the neighbours")
        s.run_tsne(channels)
        step("tSNE")
        s.run_flowsom(channels, n_clusters=3, grid=5)
        assert set(np.unique(ds.derived["FlowSOM"])) == {1, 2, 3}
        step("FlowSOM, 3 clusters")
        content = figures.MapContent(
            ds.derived["UMAP1"], ds.derived["UMAP2"], "UMAP1", "UMAP2",
            color=figures.categories_of(ds.derived["FlowSOM"], "Cluster"))
        for ext in ("svg", "pdf", "png", "tiff"):
            path = os.path.join(folder, "f." + ext)
            figures.save(content, figures.FigureStyle(), path)
            assert os.path.getsize(path) > 1000, ext
        step("figure as SVG, PDF, PNG and TIFF")

        if not a.no_window:
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
            os.environ.setdefault("MILC_SETTINGS", os.path.join(folder, "settings.ini"))
            from PyQt6 import QtWidgets

            from .ui import prepare, worker
            app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
            prepare(app)
            from .ui.shell import MainShell
            window = MainShell()
            window.open_folder(folder)
            window.samples.concatenate()
            for _ in range(3):
                worker.wait_all()
                app.processEvents()
            window.explore.clusters.setValue(3)
            window.explore.grid.setValue(5)
            window.explore.run(["UMAP", "FlowSOM"])
            for _ in range(3):
                worker.wait_all()
                app.processEvents()
            assert "FlowSOM" in s.dataset.derived or "FlowSOM" in window.session.dataset.derived
            assert window.clusters.heatmap.table is not None
            window.close()
            step("the window: folder, concatenation, UMAP and FlowSOM off the main thread")
    print(f"MIL-C {__version__}: self-test passed")
    return 0


def cmd_stats(a) -> int:
    from . import stats
    from .session import Session
    session = Session()
    session.open_folder(a.data)
    if a.wsp:
        session.set_workspace(a.wsp)
    session.set_extend_gates(not a.listed_only)
    if a.pattern:
        session.annotate(a.pattern)
    try:
        reasons = _clean(session, a)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    _say_removed(session.quality_record([session.load(e) for e in session.usable()
                                         if e.include]), reasons)
    available = session.stat_channels()
    wanted = [w.strip().lower() for w in (a.channels or "").split(",") if w.strip()]
    channels = []
    for w in wanted:
        hits = [c for c in available
                if w in (c[0].lower(), c[1].lower(), c[0].lower().split()[0],
                         c[1].lower().replace("comp-", ""))]
        if len(hits) != 1:
            print(f"{w!r}: " + ("matches " + ", ".join(h[0] for h in hits) if hits
                                else "no such parameter"), file=sys.stderr)
            return 2
        channels += hits
    by_lower = {name.lower(): name for name in stats.STATISTICS}
    chosen = []
    for w in (a.stats or "").split(","):
        if w.strip():
            if w.strip().lower() not in by_lower:
                print(f"{w!r}: choose among " + ", ".join(stats.STATISTICS), file=sys.stderr)
                return 2
            chosen.append(by_lower[w.strip().lower()])
    table = session.statistics(
        None if a.population is None
        else [_under_clean(session, p) for p in a.population.split(";")], channels, chosen)
    table.to_csv(a.out)
    print(f"{len(table.rows):,d} rows, {len(table.columns)} columns -> {a.out}")
    if a.prism:
        from . import prism
        value = a.prism_value
        if value not in table.columns:
            print(f"{value!r}: choose --prism-value among "
                  + ", ".join(table.columns[table.columns.index("Count"):]), file=sys.stderr)
            return 2
        tables = prism.from_stats_table(table, value, a.group_by.split(","), a.within)
        prism.write_pzfx(a.prism, tables)
        print(f"{len(tables)} Prism table(s) of {value} -> {a.prism}")
    return 0


def cmd_export(a) -> int:
    from .session import Session
    session = Session()
    session.open_folder(a.data)
    if a.wsp:
        session.set_workspace(a.wsp)
    session.set_extend_gates(not a.listed_only)
    try:
        reasons = _clean(session, a)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    population = _under_clean(session, a.population or "")
    entries = [e for e in session.usable() if e.include and
               (not population or population in {n.path for n in e.gates.walk()})]
    if not entries:
        print(f"no sample has the population {population!r}; there are: "
              + ", ".join(session.populations()), file=sys.stderr)
        return 2
    _say_removed(session.quality_record([session.load(e) for e in entries]), reasons)
    written = session.export_population(entries, population, a.out, a.format,
                                        compensated=not a.raw)
    print(f"{len(written)} file(s) -> {a.out}")
    return 0


# --- quality of the acquisition ------------------------------------------------

_QC_COLUMNS = ["Sample", "Events", "Checked", "Removed", "Removed %", "Intervals",
               "Time unit", "Flow rate", "Signal", "Drift", "Between", "Parameters", "Note"]


def _quality_arguments(s, apply: bool) -> None:
    """What `qc.check` can be asked, for a command: `qc` itself, or one
    that cleans the samples first (`--qc`)."""
    if apply:
        s.add_argument("--qc", action="store_true",
                       help="first remove from every sample the stretches of its "
                            "acquisition that went wrong, as “Apply to every sample” "
                            "does in the window: “Clean” is then the first population, "
                            "and the gates are under it")
    s.add_argument("--mad", type=float,
                   help="how many MADs from the acquisition's median a bin has to lie "
                        "to be suspect (6); larger removes less")
    s.add_argument("--consecutive", type=int,
                   help="good bins between two removed stretches are kept only when "
                        "there are at least this many (5)")
    s.add_argument("--events-per-bin", type=int,
                   help="by default 1/500 of the file and at least 150")
    s.add_argument("--drift", type=float,
                   help="the smallest lasting shift of a population's level that is "
                        "taken for a drift, as a share of the spread of the parameter's "
                        "events (0.25); larger removes less, 0 does not look for drifts")
    s.add_argument("--qc-channels", *(() if apply else ("--channels",)), dest="qc_channels",
                   help="the parameters to judge, comma-separated: detectors or stains; "
                        "every one but Time without it")


def _quality_settings(a) -> dict:
    given = {k: getattr(a, k, None) for k in ("mad", "consecutive", "events_per_bin", "drift")}
    return {k: v for k, v in given.items() if v is not None}


def _quality_channels(session, words: str | None) -> list[str] | None:
    """Detectors from what was typed, as `Dataset.resolve` reads it."""
    from .compensation import strip_prefix
    entries = session.usable()
    if not words or not entries:
        return None
    first, out = entries[0], []
    for w in filter(None, (w.strip() for w in words.split(","))):
        want = strip_prefix(w)[0].lower()
        hits = [c for c, m in zip(first.channels, first.markers)
                if want in (c.lower(), m.lower(), (m.lower().split() or [""])[0])]
        if len(hits) != 1:
            raise ValueError(f"{w!r}: " + ("matches " + ", ".join(hits) if hits
                                           else "no such parameter"))
        out += hits
    return out


def _check_all(session, a, keep_samples: bool = True) -> list:
    """(entry, `qc.Result`) for every file that can be read."""
    channels = _quality_channels(session, getattr(a, "qc_channels", None))
    out = []
    for e in session.usable():
        out.append((e, session.quality(e, channels=channels, **_quality_settings(a))))
        if not keep_samples:
            e.sample = None
    return out


def _time_node(session):
    from .gating import time_node
    return next((node for node in (time_node(e.gates) for e in session.usable())
                 if node is not None), None)


def _gated(session) -> list:
    """The files with gates of their own: the population of quality
    control, which every file has once it is applied, is not one — unless
    it is all there is."""
    def drawn(e):
        return [n for n in e.gates.walk() if n.gate is not None and n is not clean]
    clean = _time_node(session)
    usable = [e for e in session.usable() if e.gated]
    return [e for e in usable if drawn(e)] or usable


def _clean(session, a) -> dict[str, str] | None:
    """`--qc`: check every sample and apply what is found, through the
    session, for this run only — the project beside the files is not
    written. Where the project already has it applied, that is what is
    used. Returns why each file that could not be checked could not, or
    None where nothing was asked."""
    if not getattr(a, "qc", False):
        return None
    had = _time_node(session)
    if had is not None:
        print(f"quality control is already applied in this folder's project ({had.name!r}): "
              "it is used as it is, and nothing is checked again")
        return {}
    results = _check_all(session, a)
    judged = (next((r.channels for _, r in results if r.checked), None)
              if getattr(a, "qc_channels", None) else None)
    session.apply_quality([(e, r.intervals) for e, r in results if r.checked],
                          {**_quality_settings(a), "channels": judged}, save=False)
    return {e.name: r.reason for e, r in results if not r.checked}


def _under_clean(session, population: str | None) -> str | None:
    """A gate's path as the workspace has it, where quality control has
    put its population above every gate."""
    node = _time_node(session)
    if node is None or population is None or (population + "/").startswith(node.name + "/"):
        return population
    return "/".join(filter(None, [node.name, population]))


def _say_removed(record: dict | None, reasons: dict[str, str] | None) -> None:
    """How much each sample lost, whenever its events are the clean ones."""
    from .qc_plots import stretches
    if not record:
        return
    s = record.get("settings")
    print(f"quality control, population {record['population']!r}: "
          + (", ".join(f"{k.replace('_', ' ')} {v:g}" for k, v in s.items()
                       if k != "channels" and v is not None)
             + (f"; on {', '.join(s['channels'])}" if s.get("channels") else "")
             if s else "limits not recorded"))
    events = removed = 0
    for name, v in record["samples"].items():
        events, removed = events + v["events"], removed + v["removed"]
        if v["checked"] is False:
            said = "not checked, kept whole: " + (reasons or {}).get(name, "")
        elif v["removed"]:
            said = (f"{v['removed'] / v['events']:6.1%} removed ({v['removed']:,d} events; "
                    f"{stretches(v['intervals'], v['timestep'])})")
        else:
            said = "  0.0% removed"
        print(f"  {name:12s} {v['events']:>9,d}  {said}")
    print(f"  removed in all: {removed:,d} of {events:,d} events "
          f"({removed / max(events, 1):.1%})")


def _quality_row(name, events, checked, removed, intervals, step, result=None,
                 note: str = "") -> list:
    from . import qc
    scale = step or 1.0
    r = result if result is not None and result.checked else None
    return [name, events, "yes" if checked else "no", removed if checked else "",
            f"{100 * removed / events:.2f}" if checked and events else "",
            "; ".join(f"{lo * scale:.6g}-{hi * scale:.6g}" for lo, hi in intervals),
            "s" if step else "as stored",
            r.counts[qc.RATE] if r else "", r.counts[qc.SIGNAL] if r else "",
            r.counts[qc.DRIFT] if r else "", r.counts[qc.ISLANDS] if r else "",
            "; ".join(r.blamed) if r else "", note]


def _quality_csv(path: str, rows: list[list]) -> None:
    import csv
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding=cluster_stats.CSV_ENCODING, newline="") as f:
        out = csv.writer(f)
        out.writerow(_QC_COLUMNS)
        out.writerows(rows)


def cmd_qc(a) -> int:
    from . import qc, qc_plots
    from .session import Session
    session = Session()
    session.open_folder(a.data)
    if a.wsp:
        session.set_workspace(a.wsp)
    try:
        results = _check_all(session, a, keep_samples=False)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if not results:
        print(f"no FCS file can be read in {a.data}", file=sys.stderr)
        return 2
    steps = {e.name: qc_plots.timestep(read_header(e.path)) for e, _ in results}
    rows, events, removed = [], 0, 0
    wide = max(len(e.name) for e in session.entries)
    print(f"{'sample':{wide}s} {'events':>9s}  {'removed':>7s}  stretches, and what found them")
    for e, r in results:
        gone = int((~r.keep).sum())
        rows.append(_quality_row(e.name, len(r.keep), r.checked, gone, r.intervals,
                                 steps[e.name], r, r.reason))
        if not r.checked:
            print(f"{e.name:{wide}s} {len(r.keep):>9,d}  {'—':>7s}  not checked: {r.reason}")
            continue
        events, removed = events + len(r.keep), removed + gone
        found = [what for what in (qc.RATE, qc.SIGNAL, qc.DRIFT) if r.counts[what]]
        blamed = ", ".join(r.blamed[:3]) + ("…" if len(r.blamed) > 3 else "")
        print(f"{e.name:{wide}s} {len(r.keep):>9,d}  {r.removed:>7.1%}  "
              + (f"{qc_plots.stretches(r.intervals, steps[e.name])}  ({', '.join(found)}"
                 f"{': ' + blamed if blamed else ''})" if r.intervals else ""))
    for e in session.entries:
        if e.problem:
            rows.append(_quality_row(e.name, e.events, False, 0, [], 0.0, None, e.problem))
            print(f"{e.name:{wide}s} {e.events:>9,d}  {'—':>7s}  not checked: {e.problem}")
    checked = sum(r.checked for _, r in results)
    print(f"{checked} of {len(rows)} file(s) checked; {removed:,d} of {events:,d} events "
          f"({removed / max(events, 1):.1%}) would be removed. Nothing was: "
          "--qc on concat, stats and export removes them for that run.")
    if _time_node(session) is not None:
        print("note: this folder's project already has quality control applied; what is "
              "above is a fresh check of the whole files")
    if a.out:
        _quality_csv(a.out, rows)
        print(f"{len(rows)} rows -> {a.out}")
    if a.plots:
        try:
            written = qc_plots.save(a.plots, [(e.name, r, steps[e.name]) for e, r in results],
                                    a.plot_format, a.plot_signals)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            return 2
        print(f"{len(results)} time plot(s) -> "
              + (written[0] if len(written) == 1 else a.plots))
    return 0


def cmd_gui(a) -> int:
    # imported by name, here, where a packager can see it. Asked for through
    # __import__ it was left out of the first bundle, which then passed its
    # self-test and did not open.
    from .app import main as start
    return start(a.folder, a.quit_after)


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
        s.add_argument("--pattern",
                       help="regular expression read from each file name; every "
                            "named group is a column: '(?P<Day>\\d+)(?P<Group>[A-Za-z]+)'. "
                            "Without it there is one column, Group, guessed from "
                            "the names or as typed in the window")

    s = sub.add_parser("counts", help="every gate, against FlowJo's count")
    common(s)
    _quality_arguments(s, apply=True)
    s.set_defaults(run=cmd_counts)

    s = sub.add_parser("concat", help="balanced downsampling and concatenation")
    common(s)
    s.add_argument("--population", help="path of the gate, e.g. Cells/Live/CD45")
    s.add_argument("--mode", choices=sampling.MODES, default=sampling.EQUAL_PER_SAMPLE)
    s.add_argument("--n", type=int)
    s.add_argument("--min-events", type=int, default=1)
    s.add_argument("--group-by", default="Group",
                   help="column(s) that define a group, comma-separated")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--normalize-by", metavar="COLUMN",
                   help="batch normalisation: the column that says which batch a sample is in")
    s.add_argument("--reference", metavar="COLUMN=VALUE",
                   help="the samples that should look alike in every batch, e.g. Group=Ctrl")
    s.add_argument("--normalize-clusters", type=int, default=1,
                   help="fit within this many FlowSOM clusters; 1 aligns each marker whole")
    s.add_argument("--normalize-to", metavar="BATCH",
                   help="the batch the others are brought onto, left as acquired itself; "
                        "default: the mean of the batches")
    s.add_argument("--dry-run", action="store_true", help="the plan only")
    _quality_arguments(s, apply=True)
    s.add_argument("--out", default="concat")
    s.set_defaults(run=cmd_concat)

    s = sub.add_parser("analyze", help="UMAP, tSNE and FlowSOM on a concatenated file")
    s.add_argument("dataset", help="the .npz from concat, or any concatenated .fcs")
    s.add_argument("--channels", help="detectors or stains, comma-separated")
    s.add_argument("--transform", choices=("biex", "arcsinh", "logicle"), default="biex",
                   help="biex is what FlowJo feeds its plugins; logicle is the standard "
                        "one (T 262144, W 0.5, M 4.5)")
    s.add_argument("--cofactor", type=float, default=150.0)
    s.add_argument("--umap", action="store_true")
    s.add_argument("--neighbors", type=int, default=15)
    s.add_argument("--min-dist", type=float, default=0.5)
    s.add_argument("--tsne", action="store_true")
    s.add_argument("--perplexity", type=float, default=30.0)
    s.add_argument("--flowsom", action="store_true")
    s.add_argument("--clusters", type=int, default=5)
    s.add_argument("--grid", type=int, default=10)
    s.add_argument("--phenograph", action="store_true",
                   help="clusters without choosing how many")
    s.add_argument("--k", type=int, default=30, help="PhenoGraph's neighbours")
    s.add_argument("--resolution", type=float, default=1.0,
                   help="above 1, more and smaller clusters")
    s.add_argument("--method", choices=("leiden", "louvain"), default="leiden",
                   help="how PhenoGraph's graph is cut into communities")
    s.add_argument("--min-cluster-size", type=int, default=10,
                   help="PhenoGraph clusters of fewer events join the cluster "
                        "they are most joined to; 1 keeps them all")
    s.add_argument("--pca", action="store_true")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--group-by", default="Group")
    s.add_argument("--sample-column", default="SampleID",
                   help="for a .fcs: the column that says which sample")
    s.add_argument("--out", default="analysis")
    s.set_defaults(run=cmd_analyze)

    s = sub.add_parser("stats", help="counts, percentages and statistics of every population")
    s.add_argument("data")
    s.add_argument("--wsp", help="a FlowJo workspace; without it, the one in the folder")
    s.add_argument("--listed-only", action="store_true")
    s.add_argument("--pattern", help="columns from the file names, as for concat")
    s.add_argument("--population", help="gate paths, separated by ';'; all of them without it")
    s.add_argument("--channels", help="parameters to report on: CD8,CD4,FSC-A")
    s.add_argument("--stats", default="Median",
                   help="comma-separated: " + ", ".join(
                       ["Median", "Mean", "Geometric mean", "SD", "CV", "Robust SD",
                        "Robust CV", "MAD"]))
    s.add_argument("--out", default="statistics.csv")
    s.add_argument("--prism", metavar="FILE.pzfx",
                   help="also a GraphPad Prism file: one value per sample, by group")
    s.add_argument("--prism-value", default="% of parent", help="the column to put in it")
    s.add_argument("--group-by", default="Group", help="column(s) that make the groups")
    s.add_argument("--within", help="a table for each value of this column")
    _quality_arguments(s, apply=True)
    s.set_defaults(run=cmd_stats)

    s = sub.add_parser("export", help="the events of one population, a file per sample")
    s.add_argument("data")
    s.add_argument("--wsp", help="a FlowJo workspace; without it, the one in the folder")
    s.add_argument("--listed-only", action="store_true")
    s.add_argument("--population", help="the gate's path; every event without it")
    s.add_argument("--format", choices=("fcs", "csv"), default="fcs")
    s.add_argument("--raw", action="store_true",
                   help="as acquired, not compensated; an FCS file then carries the matrix")
    s.add_argument("--out", default="exported")
    _quality_arguments(s, apply=True)
    s.set_defaults(run=cmd_export)

    s = sub.add_parser("qc", help="the acquisition of every sample: which stretches of "
                                  "time went wrong, and how much they hold")
    s.add_argument("data")
    s.add_argument("--wsp", help="a FlowJo workspace, for its compensation; without it, "
                                 "the one in the folder")
    _quality_arguments(s, apply=False)
    s.add_argument("--out", metavar="FILE.csv", help="the table, as CSV")
    s.add_argument("--plots", metavar="PATH",
                   help="the flow rate and the signal over time, with what would be removed "
                        "shaded: FILE.pdf for one PDF with a page for each sample, anything "
                        "else is a folder that gets a file for each")
    s.add_argument("--plot-format", choices=("png", "pdf"), default="png",
                   help="of the files in that folder")
    s.add_argument("--plot-signals", type=int, default=3, metavar="N",
                   help="how many populations to draw under the flow rate, those suspect "
                        "for longest first")
    s.set_defaults(run=cmd_qc)

    s = sub.add_parser("selftest", help="the whole route on invented data: is this build sound?")
    s.add_argument("--no-window", action="store_true", help="skip the part that needs Qt")
    s.add_argument("--log", metavar="FILE",
                   help="write the steps, and any failure, to FILE; never raise")
    s.set_defaults(run=cmd_selftest)

    s = sub.add_parser("gui", help="the window (also what no command at all opens)")
    s.add_argument("folder", nargs="?")
    s.add_argument("--quit-after", type=float, metavar="SECONDS", help=argparse.SUPPRESS)
    s.set_defaults(run=cmd_gui)

    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        argv = ["gui"]
    a = ap.parse_args(argv)
    return a.run(a)


if __name__ == "__main__":
    sys.exit(main())
