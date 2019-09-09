context("Correlation matrices, possibly with a zero pattern on the Cholesky factor")

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

test_that("the size of the sample is correct", {
  N <- 10
  p <- 5
  d <- 0.25

  check_sample_size <- function(N, ...) {
    sample <- chol_mh(N = N, ...)
    expect_equal(dim(sample)[3], N)
    sample <- chol_iid(N = N, ...)
    expect_equal(dim(sample)[3], N)
  }

  # no zeros
  check_sample_size(N = N, p = p)

  # with a percentage of zeros
  check_sample_size(N = N, p = p, d = d)

  # with a predefined pattern of zeros
  dag <- rgraph(p = p, d = d, dag = TRUE)
  check_sample_size(N = N, dag = dag)
})

test_that("matrix dimension is correct", {
  N <- 10
  p <- 5
  d <- 0.25

  check_matrix_dim <- function(N, p_exp, ...) {
    sample <- chol_mh(N = N, ...)
    expect_equal(dim(sample)[1], dim(sample)[2])
    expect_equal(dim(sample)[1], p_exp)
    sample <- chol_iid(N = N, ...)
    expect_equal(dim(sample)[1], dim(sample)[2])
    expect_equal(dim(sample)[1], p_exp)
  }

  # no zeros
  check_matrix_dim(N = N, p_exp = p, p = p)

  # with a percentage of zeros
  check_matrix_dim(N = N, p_exp = p, p = p, d = d)

  # with a predefined pattern of zeros
  dag <- rgraph(p = p, d = d, dag = TRUE)
  check_matrix_dim(N = N, p_exp = p, dag = dag)
})


test_that("matrices are symmetric positive definite", {
  p <- 5
  d <- 0.25

  check_spd <- function(...) {
    sample <- chol_mh(...)
    expect_equal(sample[, , 1], t(sample[, , 1]))
    expect_gt(min(eigen(sample[, , 1])$values), 0)

    sample <- chol_iid(...)
    expect_equal(sample[, , 1], t(sample[, , 1]))
    expect_gt(min(eigen(sample[, , 1])$values), 0)
  }

  # no zeros
  check_spd(p = p)

  # with a percentage of zeros
  check_spd(p = p, d = d)

  # with a predefined zero pattern
  dag <- rgraph(p = p, d = d, dag = TRUE)
  check_spd(dag = dag)
})

test_that("matrices are of correlation", {
  p <- 5
  d <- 0.25

  check_cor <- function(...) {
    check_cor_sample <- function(sample) {
      p <- dim(sample)[1]
      for (i in 1:(p - 1)) {
        expect_equal(sample[i, i, 1], 1)
        for (j in (i + 1):p) {
          expect_gt(sample[i, j, 1], -1)
          expect_lt(sample[i, j, 1], 1)
        }
      }
    }

    sample <- chol_mh(...)
    check_cor_sample(sample)
    sample <- chol_iid(...)
    check_cor_sample(sample)
  }

  # no zeros
  check_cor(p = p)

  # with a percentage of zeros
  check_cor(p = p, d = d)

  # with a predefined zero pattern
  dag <- rgraph(p = p, d = d, dag = TRUE)
  check_cor(dag = dag)
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
