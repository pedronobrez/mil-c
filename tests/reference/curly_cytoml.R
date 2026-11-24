# CytoML on a FlowJo workspace with curly quadrants: which events it puts
# in every population, and the polygon it draws for each curly corner.
# Usage: Rscript curly_cytoml.R <workspace.wsp> <folder with the .fcs> <out folder>
#
# Writes `cytoml_members.csv` (one column per population, one row per
# event of the file, 1 where CytoML has the event in it),
# `cytoml_counts.csv` (population, CytoML's count, the count the workspace
# holds) and `cytoml_polygons.csv` (population, vertex, x, y: the polygon
# cytolib drew for each curly corner, as flowWorkspace gives it back).
suppressMessages({library(CytoML); library(flowWorkspace)})
args <- commandArgs(trailingOnly = TRUE)
ws <- open_flowjo_xml(args[1])
gs <- flowjo_to_gatingset(ws, name = 1, path = args[2])
gh <- gs[[1]]
nodes <- gs_get_pop_paths(gs, path = "full")[-1]
members <- sapply(nodes, function(p) as.integer(gh_pop_get_indices(gh, p)))
colnames(members) <- sub("^/", "", nodes)
write.csv(members, file.path(args[3], "cytoml_members.csv"), row.names = FALSE)
stats <- gh_pop_compare_stats(gh)
write.csv(data.frame(population = sub("^/", "", stats$node), cytoml = stats$openCyto.count,
                     workspace = stats$xml.count),
          file.path(args[3], "cytoml_counts.csv"), row.names = FALSE)
rows <- list()
for (p in nodes) {
  g <- gh_pop_get_gate(gh, p)
  if (is(g, "polygonGate")) {
    b <- g@boundaries
    rows[[length(rows) + 1]] <- data.frame(population = sub("^/", "", p),
                                           vertex = seq_len(nrow(b)), x = b[, 1], y = b[, 2])
  }
}
if (length(rows))
  write.csv(do.call(rbind, rows), file.path(args[3], "cytoml_polygons.csv"), row.names = FALSE)
cat("CytoML", as.character(packageVersion("CytoML")), "flowWorkspace",
    as.character(packageVersion("flowWorkspace")), "cytolib",
    as.character(packageVersion("cytolib")), "\n")
