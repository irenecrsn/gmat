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

#' Compute the anti transpose of a matrix (transpose with respect to the off-diagonal)
#'
#' @param m square matrix to compute the anti transpose
#'
#' @return The anti-transpose of m
#' @export
anti_t <- function(m) {
  p <- nrow(m)
  j <- matrix(ncol = p, nrow = p, data = 0)

  for (i in 1:p) {
    j[i, p - i + 1] <- 1
  }

  return(j %*% t(m) %*% j)
}

#' Vectorize a sample of covariance/correlation matrices
#'
#' @param sample Array, the `p x p x N` sample to vectorize
#'
#' @details Note that if the sample is of covariance matrices, as returned by [port()] and [diagdom()], the diagonal is omitted from the vectorization process.
#' @return A `p*(p - 1)/2 x N` matrix containing the vectorized sample
#' @export
vectorize <- function(sample) {
  vec_sample <- apply(sample, MARGIN = 3, function(m) {
    return(m[upper.tri(m)])
  })
  return(t(vec_sample))
}


#' Set the condition number of the matrices in a sample of covariance/correlation matrices
#'
#' @param sample Array, the `p x p x N` matrix sample
#' @param k Condition number to be set
#'
#' @return A `p x p x N` array containing the matrices with the fixed condition number
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
