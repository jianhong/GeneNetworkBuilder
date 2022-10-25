#' Subset a polished network
#' @description Subset the output of polishNetwork by a list of nodes name
#' @param graph A graphNEL object. The output of polishNetwork.
#' @param genes A list of nodes names
#' @return An object of graph.
#' @import graph
#' @export
#' @examples 
#' library(graph)
#' set.seed(123)
#' g1 <- randomEGraph(LETTERS[seq.int(15)], edges=100)
#' g1 <- subsetNetwork(g1, LETTERS[seq.int(5)])
#' plot(g1)
subsetNetwork <- function(graph, genes){
  stopifnot(
    "graph must be an object of graph"=
      inherits(graph,
               c("graphNEL", "distGraph", "clusterGraph",
                 "graphBAM", "MultiGraph")))
  stopifnot("genes must be a character vector"=is(genes, "character"))
  genes <- genes[genes %in% nodes(graph)]
  graph <- subGraph(genes, graph)
  ## remove the single nodes
  edges <- edgeWeights(graph)
  el <- lengths(edges)>0
  ep <- unlist(lapply(edges, function(.ele) names(.ele)))
  genes <- unique(c(names(edges)[el], ep))
  graph <- subGraph(genes, graph)
}
