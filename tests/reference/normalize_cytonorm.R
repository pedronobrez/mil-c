# The CytoNorm package (saeyslab/CytoNorm) on the events written by
# normalize_compare.py. Usage: Rscript normalize_cytonorm.R <work folder>
#
# The events are already on the scale they are to be normalised on, so no
# transformList is given. Each run of the manifest is either
#   clusters == 1: QuantileNorm.train / QuantileNorm.normalize, the
#                  package's own functions for normalising without clusters
#   clusters  > 1: CytoNorm.train / CytoNorm.normalize, with FlowSOM as the
#                  run says (grid, metaclusters, seed) and scale = FALSE
# with nQ = 99 and the goal the run names ("mean" or a batch).
suppressMessages({library(flowCore); library(CytoNorm); library(jsonlite)})
dir <- commandArgs(trailingOnly = TRUE)[1]
manifest <- fromJSON(file.path(dir, "manifest.json"), simplifyVector = FALSE)

read_matrix <- function(path, names) {
  v <- readBin(path, "numeric", n = file.size(path) / 4, size = 4, endian = "little")
  m <- matrix(v, ncol = length(names), byrow = TRUE)
  colnames(m) <- names
  m
}
fcs <- file.path(dir, "fcs")
unlink(fcs, recursive = TRUE)
dir.create(fcs)

for (d in manifest$datasets) {
  channels <- unlist(d$channels)
  batches <- unlist(d$batches)
  ref <- file.path(fcs, paste0(d$name, "_ref_", batches, ".fcs"))
  app <- file.path(fcs, paste0(d$name, "_x_", batches, ".fcs"))
  for (i in seq_along(batches)) {
    suppressWarnings(write.FCS(flowFrame(read_matrix(
      file.path(dir, paste0(d$name, "_ref_", batches[i], ".bin")), channels)), ref[i]))
    suppressWarnings(write.FCS(flowFrame(read_matrix(
      file.path(dir, paste0(d$name, "_x_", batches[i], ".bin")), channels)), app[i]))
  }
  for (run in d$runs) {
    out <- file.path(fcs, paste0("out_", d$name, "_", run$name))
    dir.create(out)
    files <- c(ref, app)
    labels <- c(batches, batches)
    if (run$clusters == 1) {
      model <- suppressWarnings(QuantileNorm.train(
        files = ref, labels = batches, channels = channels, transformList = NULL,
        nQ = 99, goal = run$goal, truncate_max_range = FALSE))
      suppressWarnings(QuantileNorm.normalize(
        model, files, labels, transformList = NULL, transformList.reverse = NULL,
        outputDir = out, prefix = "Norm_", truncate_max_range = FALSE))
      found <- 1
    } else {
      model <- suppressWarnings(CytoNorm.train(
        files = ref, labels = batches, channels = channels, transformList = NULL,
        outputDir = file.path(out, "train"),
        FlowSOM.params = list(nCells = 1e6, xdim = run$grid, ydim = run$grid,
                              nClus = run$clusters, scale = FALSE),
        normParams = list(nQ = 99, goal = run$goal), seed = run$seed,
        truncate_max_range = FALSE))
      suppressWarnings(CytoNorm.normalize(
        model, files, labels, transformList = NULL, transformList.reverse = NULL,
        outputDir = out, prefix = "Norm_", clean = TRUE, truncate_max_range = FALSE))
      found <- length(unique(model$fsom$metaclustering))
    }
    for (f in files) {
      ff <- read.FCS(file.path(out, paste0("Norm_", basename(f))), truncate_max_range = FALSE)
      m <- exprs(ff)[, channels, drop = FALSE]
      writeBin(as.vector(t(m)), file.path(dir, paste0(run$name, "__",
               sub("\\.fcs$", "", basename(f)), ".bin")), size = 4, endian = "little")
    }
    cat(d$name, run$name, "metaclusters", found, "\n")
  }
}
cat("CytoNorm", as.character(packageVersion("CytoNorm")), "FlowSOM",
    as.character(packageVersion("FlowSOM")), "flowCore", as.character(packageVersion("flowCore")),
    R.version.string, "\n")
