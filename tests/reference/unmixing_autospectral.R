# AutoSpectral's own functions on the same matrices as milc.unmixing.
#
#   Rscript unmixing_autospectral.R <AutoSpectral/R folder> <work folder>
#
# Reads <work>/spectra.csv (fluorochromes x detectors, first column the
# names) and <work>/raw.csv (events x detectors), and writes beside them
# what AutoSpectral gives: ols.csv, wls_cell.csv (one weighted fit per
# event), similarity.csv and complexity.txt.
#
# The functions are AutoSpectral's (github.com/DrCytometer/AutoSpectral,
# AGPL-3.0), sourced from a copy of its R folder that is not part of this
# repository: calculate_condition_number.R, cosine_similarity.R,
# unmix_ols.R, unmix_wls_per_cell.R. They need nothing but base R.

args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
work <- args[2]
for (f in c("calculate_condition_number.R", "cosine_similarity.R",
            "unmix_ols.R", "unmix_wls_per_cell.R")) {
  source(file.path(src, f))
}

spectra <- as.matrix(read.csv(file.path(work, "spectra.csv"), row.names = 1,
                              check.names = FALSE))
raw <- as.matrix(read.csv(file.path(work, "raw.csv"), check.names = FALSE))
colnames(raw) <- colnames(spectra)

write.csv(unmix.ols(raw, spectra), file.path(work, "ols.csv"), row.names = FALSE)

cells <- t(apply(raw, 1, function(y) unmix.wls.per.cell(matrix(y, nrow = 1), spectra)))
write.csv(cells, file.path(work, "wls_cell.csv"), row.names = FALSE)

write.csv(cosine.similarity(spectra), file.path(work, "similarity.csv"))
# rounded to 2 places by AutoSpectral; the unrounded number beside it
sv <- svd(spectra)$d
writeLines(c(format(calculate.condition.number(spectra)),
             format(max(sv) / min(sv), digits = 15)),
           file.path(work, "complexity.txt"))
cat("AutoSpectral functions ran on", nrow(raw), "events,", nrow(spectra), "spectra\n")
