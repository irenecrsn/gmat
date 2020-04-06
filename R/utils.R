#' Random generation of acyclic digraphs and undirected graphs
#'
#' Wrapper of functionality from package `igraph` for random generation of
#' graphs.
#'
#' @rdname rgraph
#'
#' @param p Number of vertices of the sampled graph.
#' @param d Proportion of edges in the generated graph.
#' @param dag Whether the generated graph should be acyclic directed.
#' @param ordered When generating an acyclic directed graph, whether the nodes
#'  should follow the ancestral order `1, ..., p`.
#'
#' @details When `dag = FALSE`, the graph is sampled from an Erdos-Renyi model.
#' In the case where `dag = TRUE`, the upper triangle of the adjacency matrix of
#' an Erdos-Renyi model is taken as the adjacency matrix for the acyclic
#' digraph. This preserves the proportion of edges `d`.
#'
#' @return The generated graph.
#'
#' @examples
#' ## Random undirected graph with 3 nodes and 50% density of edges
#' rgraph(p = 3, d = 0.5)
#'
#' ## Random directed acyclic graphs
#' # Following the natural ancestral order 1, ..., p
#' dag <- rgraph(p = 6, d = 0.5, dag = TRUE)
#' igraph::topo_sort(dag)
#'
#' # Following a random ancestral order
#' dag <- rgraph(p = 6, d = 0.5, dag = TRUE, ordered = FALSE)
#' igraph::topo_sort(dag)
#' @export
rgraph <- function(p, d, dag = FALSE, ordered = TRUE) {
  g <- igraph::sample_gnp(n = p, p = d)

  if (dag == TRUE) {
    dag_am <- igraph::as_adjacency_matrix(g, type = "upper")
    if (ordered == FALSE) {
      random_order <- sample(seq.int(from = 1, to = p))
      dag_am <- dag_am[random_order, random_order]
    }
    g <- igraph::graph_from_adjacency_matrix(dag_am, mode = "directed")
  }

  return(g)
}


#' Get the upper factor of the upper Cholesky decomposition of a symmetric
#' positive definite matrix.
#'
#' @param m Matrix to factorize.
#'
#' @details The upper factor `U` such that `m = U %*% t(U)`. `U` is equal to the
#' transpose with respect to the anti-diagonal of the standard Cholesky factor
#' `L` in `m_rev = L %*% t(L)`, where `m_rev` is the matrix resulting from
#' reverting the order of rows and columns in `m` (see Córdoba et al., 2019,
#' Section 2.2 for more details). The function uses the base `chol` method.
#' @return A `p*p` upper triangular matrix.
#'
#' @references Córdoba I., Varando G., Bielza C., Larrañaga P., Generating
#' random Gaussian graphical models, _arXiv_:1909.01062, 2019.
#' @export
uchol <- function(m) {
  p <- nrow(m)

  # Reverse order of rows and columns of the given matrix
  m <- m[p:1, p:1]

  # Upper standard Cholesky factor
  u <- chol(m)

  # Transpose with respect to the antidiagonal
  return(t(u[p:1, p:1]))
}

#' Vectorize a sample of covariance/correlation matrices
#'
#' @param sample Array, the `p x p x N` sample to vectorize.
#'
#' @details Note that if the sample is of covariance matrices, as returned by [port()] and [diagdom()], the diagonal is omitted from the vectorization process.
#' @return A `p*(p - 1)/2 x N` matrix containing the vectorized sample.
#' @export
vectorize <- function(sample) {
  vec_sample <- apply(sample, MARGIN = 3, function(m) {
    return(m[upper.tri(m)])
  })
  return(t(vec_sample))
}


#' Set the condition number of the matrices in a sample of covariance/correlation matrices
#'
#' @param sample Array, the `p x p x N` matrix sample.
#' @param k Condition number to be set.
#'
#' @return A `p x p x N` array containing the matrices with the fixed condition
#' number.
#' @export
set_cond_number <- function(sample, k) {
  N <- dim(sample)[3]
  p <- dim(sample)[1]

  for (n in 1:N) {
    eig_val <- eigen(sample[, , n])$values
    delta <- (max(eig_val) - k * min(eig_val)) / (k - 1)
    sample[, , n] <- sample[, , n] + diag(x = delta, nrow = p)
  }

  return(sample)
}


#' Moral DAG from non chordal UG
#'
#' Find the DAG with no v-structures whose skeleton is a
#' triangulation of a given undirected graph.
#'
#' @param ug An igraph undirected graph.
#' @return An igraph acyclic directed graph orientation.
#' @export
ug_to_dag <- function(ug) {

  # We triangulate the undirected graph if it is not chordal
  ug_cover <- igraph::is_chordal(ug, newgraph = TRUE)$newgraph

  # We get the max_cardinality sort == perfect ordering
  # The perfect ordering will be the ancestral ordering of orientation
  # By construction this ordering cannot induce v-structures
  ug_cover_mcsort <- igraph::max_cardinality(ug_cover)
  dag_topo_sort <- ug_cover_mcsort$alpha
  inv <- ug_cover_mcsort$alpham1
  dag_mat <- igraph::as_adjacency_matrix(ug_cover, sparse = FALSE)
  dag_mat <- dag_mat[dag_topo_sort, dag_topo_sort]
  dag_mat[lower.tri(dag_mat)] <- 0

  dag <- igraph::graph_from_adjacency_matrix(dag_mat[inv, inv], mode = "directed")
  return(dag)
}
