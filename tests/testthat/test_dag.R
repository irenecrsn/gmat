context("Possible zero pattern on the Cholesky factor represented by a DAG")

test_that("upper Cholesky factors are of Cholesky", {
  p <- 5
  d <- 0.2

  expect_cholesky <- function(u) {
    m <- tcrossprod(u)
    expect_equal(u, uchol(m))
  }

  # No zero pattern
  u_sample <- mh_u(p = p)[, , 1]
  expect_cholesky(u_sample)

  # Random zero pattern corresponding to a "natural" dag
  dag <- rgraph(p = p, d = d, dag = TRUE)
  u_sample <- mh_u(dag = dag)[, , 1]
  expect_cholesky(u_sample)

  # Random zero pattern corresponding to an arbitrary dag
  dag <- rgraph(p = p, d = d, dag = TRUE, ordered = FALSE)
  u_sample <- mh_u(dag = dag)[, , 1]
  topsort <- gRbase::topoSort(igraph::as_graphnel(dag), index = TRUE)
  expect_cholesky(u_sample[topsort, topsort])
})


test_that("the dag structure is preserved", {
  p <- 10
  d <- 0.25

  expect_equal_dag <- function(m, dag) {
    topsort <- gRbase::topoSort(igraph::as_graphnel(dag), index = TRUE)

    madj <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
    madj <- madj[topsort, topsort]

    u <- uchol(m[topsort, topsort])
    u <- (zapsmall(u) != 0) # For ignoring numeric errors
    diag(u) <- FALSE

    expect_equal(length(which((madj - u) != 0)), 0)
  }

  # With the natural ancestral order
  dag <- rgraph(p = p, d = d, dag = TRUE)

  sample <- chol_mh(dag = dag)
  expect_equal_dag(m = sample[, , 1], dag = dag)

  sample <- chol_iid(dag = dag)
  expect_equal_dag(m = sample[, , 1], dag = dag)

  # With a random ancestral order
  dag <- rgraph(p = p, d = d, dag = TRUE, ordered = FALSE)

  sample <- chol_mh(dag = dag)
  expect_equal_dag(m = sample[, , 1], dag = dag)

  sample <- chol_iid(dag = dag)
  expect_equal_dag(m = sample[, , 1], dag = dag)
})
