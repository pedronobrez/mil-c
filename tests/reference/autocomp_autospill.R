# AutoSpill (carlosproca/autospill) and flowStats::spillover on the
# single-stain controls written by autocomp_compare.py.
# Usage: Rscript autocomp_autospill.R <work folder>
#
# For each data set of the manifest, in its folder:
#   own/     AutoSpill as its batch script `calculate_compensation_minimal.r`
#            runs it: its tessellation gate on scatter, get.marker.spillover,
#            refine.spillover. The gate is written out (gate_<dye>.csv, the
#            row numbers kept) so the Python side can use the same events.
#   all/     the same two steps on every event, no gate: the algorithms
#            alone, on exactly the events MIL-C is given
#   posnegpop_all.csv  get.marker.spillover.posnegpop on every event
#   flowstats.csv      flowStats::spillover(method = "median"), the
#            unstained control as the baseline, stains matched in order
# Spillover matrices are written as rows = dye, columns = detectors.
suppressMessages({library(flowCore); library(autospill); library(jsonlite)})
dir <- commandArgs(trailingOnly = TRUE)[1]
manifest <- fromJSON(file.path(dir, "manifest.json"), simplifyVector = FALSE)

write_matrix <- function(m, path) write.csv(m, path)

for (d in manifest$datasets) {
  out <- file.path(dir, d$name)
  stained <- Filter(function(f) !is.null(f$dye), d$files)
  ctl <- data.frame(filename = sapply(stained, `[[`, "file"),
                    dye = sapply(stained, `[[`, "dye"),
                    antigen = sapply(stained, `[[`, "dye"),
                    wavelength = NA, stringsAsFactors = FALSE)
  ctl_file <- file.path(out, "fcs_control.csv")
  write.csv(ctl, ctl_file, row.names = FALSE, na = "")

  asp <- get.autospill.param("minimal")
  asp$verbose <- FALSE
  asp$default.scatter.parameter <- unlist(d$scatter)
  asp$gate.parameter.file.name <- NULL
  asp$scatter.parameter.file.name <- NULL
  asp$transformation.parameter.file.name <- NULL
  fc <- read.flow.control(file.path(out, "fcs"), ctl_file, asp)

  runs <- list()
  t0 <- proc.time()[["elapsed"]]
  gate <- tryCatch(gate.flow.data(fc, asp), error = function(e) {
    cat("AutoSpill's gate failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(gate)) {
    runs$own <- gate
    for (s in fc$sample)
      write.csv(data.frame(row = gate[[s]]), file.path(out, paste0("gate_", s, ".csv")),
                row.names = FALSE)
  }
  runs$all <- lapply(fc$sample, function(s) seq_len(sum(fc$event.sample == s)))
  names(runs$all) <- fc$sample
  cat(sprintf("%s: gate %.1f s\n", d$name, proc.time()[["elapsed"]] - t0))

  for (r in names(runs)) {
    dest <- file.path(out, r)
    dir.create(dest, showWarnings = FALSE)
    t0 <- proc.time()[["elapsed"]]
    first <- get.marker.spillover(TRUE, runs[[r]], fc, asp)
    ref <- tryCatch(refine.spillover(first, NULL, runs[[r]], fc, asp),
                    error = function(e) { cat("refine:", conditionMessage(e), "\n"); NULL })
    took <- proc.time()[["elapsed"]] - t0
    m0 <- first$coef
    dimnames(m0) <- list(fc$marker.original, fc$marker.original)
    write_matrix(m0, file.path(dest, "initial.csv"))
    if (!is.null(ref)) {
      m <- ref$spillover
      dimnames(m) <- list(fc$marker.original, fc$marker.original)
      write_matrix(m, file.path(dest, "final.csv"))
      write.csv(ref$convergence, file.path(dest, "convergence.csv"), row.names = FALSE)
    }
    writeLines(sprintf("%.3f", took), file.path(dest, "seconds.txt"))
    cat(sprintf("%s/%s: %.1f s, %s\n", d$name, r, took,
                if (is.null(ref)) "no refinement" else
                  sprintf("%d iterations", nrow(ref$convergence))))
    if (r == "all") {
      pn <- get.marker.spillover.posnegpop(runs[[r]], fc, asp)
      dimnames(pn) <- list(fc$marker.original, fc$marker.original)
      write_matrix(pn, file.path(dest, "..", "posnegpop_all.csv"))
    }
  }

  unst <- Filter(function(f) is.null(f$dye), d$files)
  if (length(unst) == 1 && requireNamespace("flowStats", quietly = TRUE)) {
    suppressMessages(library(flowStats))
    files <- c(unst[[1]]$file, sapply(stained, `[[`, "file"))
    fs <- read.flowSet(file.path(out, "fcs", files), transformation = FALSE,
                       truncate_max_range = FALSE)
    keep <- c(unlist(d$scatter), sapply(stained, `[[`, "dye"))
    fs <- fs[, keep]
    sampleNames(fs) <- c("unstained", sapply(stained, `[[`, "dye"))
    sp <- suppressMessages(spillover(fs, unstained = "unstained", fsc = d$scatter[[1]],
                                     ssc = d$scatter[[2]], method = "median",
                                     stain_match = "regexpr", exact_match = TRUE))
    write_matrix(sp, file.path(out, "flowstats.csv"))
  }
}
cat("versions:", as.character(packageVersion("autospill")), as.character(packageVersion("flowCore")),
    if (requireNamespace("flowStats", quietly = TRUE)) as.character(packageVersion("flowStats")) else "",
    as.character(packageVersion("MASS")), "\n")
