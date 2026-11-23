# flowFit on the data it ships with: Quah & Parish's CD4 T cells labelled
# with CFSE, CellTrace Violet and Cell Proliferation Dye (their Fig. 2a).
# Usage: Rscript proliferation_flowfit.R <work folder>
#
# As the package's vignette does it: the undivided peak is fitted on the
# unstimulated sample (`parentFitting`), then the stimulated one with that
# position and width as the start (`proliferationFitting`), position,
# width and distance between generations all free. Writes the events of
# each sample (`quah_<k>.csv`: the fluorescence is log10 already) and what
# flowFit found (`flowfit.csv`).
suppressMessages({library(flowCore); library(flowFitExampleData); library(flowFit)})
dir <- commandArgs(trailingOnly = TRUE)[1]
data(QuahAndParish)
for (k in seq_along(QuahAndParish))
  write.csv(exprs(QuahAndParish[[k]]), file.path(dir, paste0("quah_", k, ".csv")),
            row.names = FALSE)

# the stimulated sample labelled with each dye, and the channel it is read on
stained <- list(c("<FITC-A>", 2), c("<APC-A>", 3), c("<Alexa Fluor 405-A>", 4))
rows <- list()
for (pair in stained) {
  channel <- pair[1]
  k <- as.integer(pair[2])
  parent <- parentFitting(QuahAndParish[[1]], channel)
  fit <- proliferationFitting(QuahAndParish[[k]], channel, parent@parentPeakPosition,
                              parent@parentPeakSize)
  percent <- unlist(getGenerations(fit))
  rows[[length(rows) + 1]] <- data.frame(
    channel = channel, sample = k, parent_position = parent@parentPeakPosition,
    parent_size = parent@parentPeakSize, position = fit@parentPeakPosition,
    size = fit@parentPeakSize, distance = fit@generationsDistance,
    index = proliferationIndex(fit), percent = paste(percent, collapse = " "))
}
write.csv(do.call(rbind, rows), file.path(dir, "flowfit.csv"), row.names = FALSE)
cat("flowFit", as.character(packageVersion("flowFit")), "flowCore",
    as.character(packageVersion("flowCore")), "\n")
