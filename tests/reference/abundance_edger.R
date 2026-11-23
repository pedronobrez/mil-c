# edgeR (and diffcyt's DA-edgeR) on the count tables written by
# abundance_compare.py. Usage: Rscript abundance_edger.R <work folder>
#
# For every table, five ways of asking the same question:
#   ql_matched   what milc/abundance.py says it imitates: one common
#                dispersion (estimateGLMCommonDisp), glmQLFit as before
#                edgeR 4 (legacy = TRUE), no abundance trend, glmQLFTest
#   ql_default   what a user of edgeR 4 gets: estimateDisp, glmQLFit,
#                glmQLFTest, all at their defaults
#   lrt_matched  estimateGLMCommonDisp + estimateGLMTagwiseDisp(prior.df =
#                10, trend = FALSE), glmFit, glmLRT
#   lrt_disp     estimateDisp(trend.method = "none"), glmFit, glmLRT: the
#                calls diffcyt::testDA_edgeR makes, written out so that the
#                dispersions can be kept
#   diffcyt      diffcyt::testDA_edgeR itself (no filter of rare clusters,
#                no normalisation factors: its default)
suppressMessages({library(edgeR); library(jsonlite)})
dir <- commandArgs(trailingOnly = TRUE)[1]
tables <- fromJSON(file.path(dir, "tables.json"), simplifyVector = FALSE)
have_diffcyt <- suppressWarnings(suppressMessages(requireNamespace("diffcyt", quietly = TRUE)))
if (have_diffcyt) suppressMessages(library(SummarizedExperiment))

rows <- list()
keep <- function(id, variant, tab, dispersion, statistic, common = NA, prior_df = NA) {
  rows[[length(rows) + 1]] <<- data.frame(
    table = id, variant = variant, cluster = seq_len(nrow(tab)) - 1,
    log2fc = tab$logFC, dispersion = dispersion, statistic = statistic,
    p = tab$PValue, fdr = p.adjust(tab$PValue, "BH"),
    common = common, prior_df = prior_df)
}

for (t in tables) {
  counts <- matrix(unlist(t$counts), nrow = length(t$counts), byrow = TRUE)
  rownames(counts) <- seq_len(nrow(counts))
  groups <- unlist(t$groups)
  design <- cbind(1, as.numeric(groups == t$contrast[[1]]))
  if (!is.null(t$covariate)) {
    cov <- factor(unlist(t$covariate), levels = unique(unlist(t$covariate)))
    if (nlevels(cov) > 1) design <- cbind(design, model.matrix(~cov)[, -1, drop = FALSE])
  }
  y <- DGEList(counts)        # lib.size = column sums, norm.factors = 1

  common <- estimateGLMCommonDisp(y, design)$common.dispersion
  fit <- glmQLFit(y, design, dispersion = common, legacy = TRUE,
                  abundance.trend = FALSE, robust = FALSE)
  keep(t$id, "ql_matched", glmQLFTest(fit, coef = 2)$table, rep(common, nrow(counts)),
       glmQLFTest(fit, coef = 2)$table$F, common, fit$df.prior[1])

  yd <- estimateDisp(y, design)
  fit <- glmQLFit(yd, design)
  res <- glmQLFTest(fit, coef = 2)
  keep(t$id, "ql_default", res$table, rep_len(fit$dispersion, nrow(counts)), res$table$F,
       yd$common.dispersion, fit$df.prior[1])

  yt <- estimateGLMTagwiseDisp(estimateGLMCommonDisp(y, design), design,
                               prior.df = 10, trend = FALSE)
  res <- glmLRT(glmFit(yt, design), coef = 2)
  keep(t$id, "lrt_matched", res$table, yt$tagwise.dispersion, res$table$LR, common, 10)

  # edgeR 4.10.7: estimateDisp(trend.method = "none") stops ("non-conformable
  # arrays") when a row has fewer than min.row.sum = 5 counts in all, so
  # those rows are left out here and come back as NA. diffcyt's own filter
  # (3 cells in half the samples) would have dropped them too.
  sub <- rowSums(counts) >= 5
  fill <- function(v) { out <- rep(NA_real_, nrow(counts)); out[sub] <- v; out }
  ys <- DGEList(counts[sub, , drop = FALSE], lib.size = colSums(counts))
  yn <- estimateDisp(ys, design, trend.method = "none")
  res <- glmLRT(glmFit(yn, design), coef = 2)
  keep(t$id, "lrt_disp",
       data.frame(logFC = fill(res$table$logFC), PValue = fill(res$table$PValue)),
       fill(yn$tagwise.dispersion), fill(res$table$LR), yn$common.dispersion, yn$prior.df[1])

  if (have_diffcyt) {
    d <- SummarizedExperiment(
      assays = list(counts = counts[sub, , drop = FALSE]),
      rowData = data.frame(cluster_id = factor(rownames(counts)[sub],
                                               levels = rownames(counts)[sub])))
    contrast <- c(0, 1, rep(0, ncol(design) - 2))
    out <- as.data.frame(rowData(diffcyt::testDA_edgeR(d, design, contrast,
                                                       min_cells = 0, min_samples = 0)))
    keep(t$id, "diffcyt", data.frame(logFC = fill(out$logFC), PValue = fill(out$p_val)),
         NA, fill(out$LR))
  }
}
write.csv(do.call(rbind, rows), file.path(dir, "edger.csv"), row.names = FALSE)
cat("edgeR", as.character(packageVersion("edgeR")), "limma", as.character(packageVersion("limma")),
    if (have_diffcyt) paste("diffcyt", as.character(packageVersion("diffcyt"))) else "diffcyt absent",
    R.version.string, "\n")
