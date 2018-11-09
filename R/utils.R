#' Random generation of acyclic digraphs and undirected graphs
#'
#' Wrapper of functionality from package `igraph` for random generation of
#' graphs. 
#' 
#' @rdname rgraph
#'
#' @param p Number of vertices of the sampled graph
#' @param d Proportion of edges in the generated graph
#' @param dag Whether the generated graph should be acyclic directed
#'
#' @details When `dag = FALSE`, the graph is sampled from an Erdos-Renyi model.
#' In the case where `dag = TRUE`, the upper triangle of the adjacency matrix of 
#' an Erdos-Renyi model is taken as the adjacency matrix for the acyclic
#' digraph. This preserves the proportion of edges `d`.
#'
#' @return g The generated graph. If `dag = TRUE`, the nodes follow the
#' ancestral order `1, ..., p`, where `1` has no parents.
#' @export
rgraph <- function(p, d, dag = FALSE) {
	
	g <- igraph::sample_gnp(n = p, p = d)

	if (dag == TRUE) {
		dag_am <- igraph::as_adjacency_matrix(g, type = "upper")
		g <- igraph::graph_from_adjacency_matrix(dag_am, mode = "directed")
	} 

	return(g)
}
