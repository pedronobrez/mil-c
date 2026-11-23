# PeacoQC and flowAI on the acquisitions written by qc_compare.py.
# Usage: Rscript qc_peacoqc.R <work folder>
#
# PeacoQC is given exactly the values milc.qc.check judged (compensated,
# each parameter on its display scale), in the file's order, at its
# defaults (MAD = 6, IT_limit = 0.6, consecutive_bins = 5): once whole
# ("all": isolation tree, then MAD) and once with the MAD step alone.
# flowAI is given the values as acquired, as it expects, and asked for its
# flow-rate and signal checks (remove_from = "FR_FS"); its third check,
# of the dynamic range, has no counterpart in ours and is left out.
suppressMessages({library(flowCore); library(PeacoQC); library(jsonlite)})
dir <- commandArgs(trailingOnly = TRUE)[1]
items <- fromJSON(file.path(dir, "manifest.json"), simplifyVector = FALSE)
have_flowai <- suppressWarnings(suppressMessages(requireNamespace("flowAI", quietly = TRUE)))

read_matrix <- function(path, n, names) {
  m <- matrix(readBin(path, "numeric", n = n * length(names), size = 4, endian = "little"),
              nrow = n, byrow = TRUE)
  colnames(m) <- names
  m
}
frame <- function(m, name, timestep = NULL) {
  ff <- flowFrame(m)
  keyword(ff)[["$FIL"]] <- name
  if (!is.null(timestep)) keyword(ff)[["$TIMESTEP"]] <- as.character(timestep)
  ff
}
write_flags <- function(good, path) writeBin(as.raw(as.integer(good)), path)

notes <- list()
for (it in items) {
  channels <- unlist(it$channels)
  time <- readBin(file.path(dir, it$time), "numeric", n = it$n, size = 8, endian = "little")
  x <- cbind(read_matrix(file.path(dir, it$scaled), it$n, channels), Time = time)
  ff <- frame(x, it$name)
  note <- ""
  for (how in c("all", "MAD")) {
    good <- tryCatch({
      res <- suppressWarnings(suppressMessages(
        PeacoQC(ff, channels, determine_good_cells = how, plot = FALSE, save_fcs = FALSE,
                report = FALSE, output_directory = NULL)))
      if (how == "all") note <- paste0(
        "bins of ", res$EventsPerBin, "; IT ", round(res$ITPercentage, 2),
        "%, MAD ", round(res$MADPercentage, 2), "%, consecutive ",
        round(res$ConsecutiveCellsPercentage, 2), "%",
        if (length(res$WeirdChannels$Changing_channel) &&
            res$WeirdChannels$Changing_channel != "No increasing or decreasing effect")
          paste0("; ", res$WeirdChannels$Changing_channel) else "")
      res$GoodCells
    }, error = function(e) { note <<- paste("PeacoQC", how, "failed:", conditionMessage(e)); NULL })
    if (!is.null(good)) write_flags(good, file.path(dir, paste0(it$name, ".peacoqc_", how, ".bin")))
  }
  if (have_flowai) {
    raw <- cbind(read_matrix(file.path(dir, it$raw), it$n, channels), Time = time)
    good <- tryCatch({
      out <- capture.output(bad <- suppressWarnings(suppressMessages(
        flowAI::flow_auto_qc(frame(raw, it$name, it$timestep), remove_from = "FR_FS",
                             output = 3, timeCh = "Time", timestep = it$timestep,
                             html_report = FALSE, mini_report = FALSE, fcs_QC = FALSE,
                             folder_results = FALSE))))
      g <- rep(TRUE, it$n); g[unlist(bad)] <- FALSE; g
    }, error = function(e) { note <<- paste(note, "| flowAI failed:", conditionMessage(e)); NULL })
    if (!is.null(good)) write_flags(good, file.path(dir, paste0(it$name, ".flowai.bin")))
  }
  notes[[it$name]] <- note
  cat(it$name, note, "\n")
}
write_json(notes, file.path(dir, "notes.json"), auto_unbox = TRUE)
cat("PeacoQC", as.character(packageVersion("PeacoQC")),
    if (have_flowai) paste("flowAI", as.character(packageVersion("flowAI"))) else "flowAI absent",
    "flowCore", as.character(packageVersion("flowCore")), R.version.string, "\n")
