# Comparisons with the original packages

Four methods of MIL-C were written from the papers. The scripts here run
each beside the software it takes after, on the same input, and print what
was measured. The results of the run of 2026-11-23 are written up in
`docs/comparacao.md`.

Nothing here is part of the test suite: no file is named `test_*.py`, and
`python3 -m pytest` does not look in this folder for tests. R and the
packages below are needed only to run these scripts.

| script | ours | beside |
|---|---|---|
| `abundance_compare.py` + `abundance_edger.R` | `milc.abundance.differential` | edgeR (`glmQLFit`/`glmQLFTest`, `glmFit`/`glmLRT`), `diffcyt::testDA_edgeR` |
| `phenograph_compare.py` | `milc.graphcluster.phenograph` | the authors' `PhenoGraph` (PyPI), with its Louvain and with Leiden |
| `qc_compare.py` + `qc_peacoqc.R` | `milc.qc.check` | PeacoQC, flowAI |
| `normalize_compare.py` + `normalize_cytonorm.R` | `milc.normalize.fit` / `apply` | CytoNorm (`saeyslab/CytoNorm`) |
| `curly_compare.py` + `curly_cytoml.R` | `milc.gating.CurlyQuadrantGate` read by `milc.wsp` | CytoML (`flowjo_to_gatingset`), event by event |

## Installing

Everything goes into a folder of your own, not into the system library.

```sh
export REF=/some/scratch/folder
export R_LIBS_USER=$REF/Rlib
mkdir -p $R_LIBS_USER

# R packages: CRAN and Bioconductor
Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); options(repos = "https://cloud.r-project.org")
  install.packages(c("BiocManager", "jsonlite"), lib = lib)
  BiocManager::install(c("edgeR", "diffcyt", "flowCore", "PeacoQC", "flowAI", "FlowSOM"),
                       lib = lib, ask = FALSE, update = FALSE)
  install.packages(c("emdist", "dplyr", "stringr", "pheatmap", "ggplot2", "ggpubr",
                     "gridExtra", "ggridges", "tidyr"), lib = lib)'

# CytoNorm, from its authors' repository
curl -sL https://github.com/saeyslab/CytoNorm/archive/refs/heads/master.tar.gz | tar xz -C $REF
R CMD INSTALL --no-docs --no-build-vignettes -l $R_LIBS_USER $REF/CytoNorm-master

# PhenoGraph, in a virtual environment that also sees the project's own dependencies
python3 -m venv --system-site-packages $REF/venv
$REF/venv/bin/pip install PhenoGraph
```

On an Apple-silicon Mac, PhenoGraph's own Louvain is an Intel program and
needs Rosetta. diffcyt and flowAI are optional: the scripts leave their
columns out when they are not installed.

Versions of the run that is written up: R 4.6.1, Bioconductor 3.23, edgeR
4.10.7, limma 3.68.5, diffcyt 1.32.1, flowCore 2.24.0, PeacoQC 1.22.0,
flowAI 1.42.0, FlowSOM 2.20.0, CytoNorm 2.0.13 (commit 515a4e9),
PhenoGraph 1.5.7 with leidenalg 0.12.0, Python 3.14.6.

## Running

From the repository's folder:

```sh
export MILC_REF_WORK=$REF/work          # where everything is written
export MILC_REF_DATA=/path/to/the/fcs   # optional: the experiment's .fcs files and .wsp
export RSCRIPT=Rscript                  # optional: which Rscript

python3 tests/reference/abundance_compare.py 100     # tables per setting
$REF/venv/bin/python tests/reference/phenograph_compare.py
python3 tests/reference/qc_compare.py
python3 tests/reference/normalize_compare.py
python3 tests/reference/curly_compare.py             # needs CytoML: BiocManager::install("CytoML")
```

On the day of the run, on a machine busy with other work, they took about
10 minutes, an hour, 50 minutes and 10 minutes: PhenoGraph's own Jaccard
step and PeacoQC on the larger files are what takes long.

Each prints its tables and leaves them in `$MILC_REF_WORK/<name>/report.md`
(and `.json`). Without `MILC_REF_DATA` only the simulated data are used.
The real data are only read; nothing is written beside them.

`MILC_REF_DATA` is expected to hold the experiment these were written
around (files named `<day><Case|Ctrl><animal>.fcs` and the FlowJo
workspace with the gates down to `Live/CD45`); `common.py` says how it is
read.

`phenograph_compare.py` calls nothing but
`graphcluster.phenograph(X, k=, resolution=, seed=)`, so it can be run
again after that function changes inside.
