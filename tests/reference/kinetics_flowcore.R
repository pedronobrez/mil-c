# Kinetics written again in R from the definitions of milc.kinetics, the
# files read by flowCore: 100 slices of Time x $TIMESTEP over the
# acquisition, the median of FITC-A and of FITC-A / PI-A in each (slices of
# at least 5 events), a centred moving average of 5 that skips empty slices;
# baseline = mean of the raw slices in the first sixth, peak = highest
# smoothed slice after it, AUC = trapezia of (smoothed - baseline).
# usage: Rscript kinetics_flowcore.R <list of files> <out.csv>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(flowCore))
files <- readLines(args[1])
out <- NULL
for (path in files) {
  ff <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  e <- exprs(ff)
  step <- suppressWarnings(as.numeric(keyword(ff, "$TIMESTEP")[[1]]))
  if (length(step) == 0 || is.na(step)) step <- 1
  t <- e[, "Time"] * step
  for (what in c("FITC-A", "ratio")) {
    v <- if (what == "ratio") e[, "FITC-A"] / e[, "PI-A"] else e[, "FITC-A"]
    keep <- is.finite(t) & is.finite(v)
    if (what == "ratio") keep <- keep & e[, "PI-A"] > 0
    tt <- t[keep]; vv <- v[keep]
    first <- min(tt); last <- max(tt); width <- (last - first) / 100
    k <- pmin(pmax(floor((tt - first) / width), 0), 99)
    y <- sapply(0:99, function(i) { s <- vv[k == i]; if (length(s) >= 5) median(s) else NA })
    centre <- first + (0:99 + 0.5) * width
    sm <- sapply(1:100, function(i) {
      if (is.na(y[i])) return(NA)
      mean(y[max(1, i - 2):min(100, i + 2)], na.rm = TRUE)
    })
    cut <- first + (last - first) / 6
    base <- mean(y[centre >= first & centre <= cut], na.rm = TRUE)
    after <- which(centre >= cut & centre <= last & !is.na(y))
    top <- after[which.max(sm[after])]
    d <- sm[after] - base
    auc <- sum(diff(centre[after]) * (head(d, -1) + tail(d, -1)) / 2)
    out <- rbind(out, data.frame(file = path, what = what, baseline = base, peak = sm[top],
                                 time_of_peak = centre[top], auc = auc,
                                 median_slice_50 = y[50]))
  }
}
write.csv(out, args[2], row.names = FALSE)
cat(nrow(out) / 2, "files\n")
