context("SPD matrices, possibly with a zero pattern")

test_that("the size of the sample is correct", {
  N <- 10
  p <- 5
  d <- 0.25

  check_sample_size <- function(N, ...) {
    sample <- port(N = N, ...)
    expect_equal(dim(sample)[3], N)

    sample <- diagdom(N = N, ...)
    expect_equal(dim(sample)[3], N)

    sample <- port_chol(N = N, ...)
    expect_equal(dim(sample)[3], N)
  }

  # no zeros
  check_sample_size(N = N, p = p)

  # with a percentage of zeros
  check_sample_size(N = N, p = p, d = d)

  # with a predefined pattern of zeros
  ug <- rgraph(p = p, d = d)
  check_sample_size(N = N, ug = ug)
})

test_that("matrix dimension is correct", {
  N <- 10
  p <- 5
  d <- 0.25

  check_matrix_dim <- function(N, p_exp, ...) {
    sample <- port(N = N, ...)
    expect_equal(dim(sample)[1], dim(sample)[2])
    expect_equal(dim(sample)[1], p_exp)

    sample <- diagdom(N = N, ...)
    expect_equal(dim(sample)[1], dim(sample)[2])
    expect_equal(dim(sample)[1], p_exp)

    sample <- port_chol(N = N, ...)
    expect_equal(dim(sample)[1], dim(sample)[2])
    expect_equal(dim(sample)[1], p_exp)
  }

  # no zeros
  check_matrix_dim(N = N, p_exp = p, p = p)

  # with a percentage of zeros
  check_matrix_dim(N = N, p_exp = p, p = p, d = d)

  # with a predefined pattern of zeros
  ug <- rgraph(p = p, d = d)
  check_matrix_dim(N = N, p_exp = p, ug = ug)
})


test_that("matrices are symmetric positive definite", {
  p <- 5
  d <- 0.25

  check_spd <- function(...) {
    sample <- port(...)
    expect_equal(sample[, , 1], t(sample[, , 1]))
    expect_gt(min(eigen(sample[, , 1])$values), 0)

    sample <- diagdom(...)
    expect_equal(sample[, , 1], t(sample[, , 1]))
    expect_gt(min(eigen(sample[, , 1])$values), 0)

    sample <- port_chol(...)
    expect_equal(sample[, , 1], t(sample[, , 1]))
    expect_gt(min(eigen(sample[, , 1])$values), 0)
  }

  # no zeros
  check_spd(p = p)

  # with a percentage of zeros
  check_spd(p = p, d = d)

  # with a predefined zero pattern
  ug <- rgraph(p = p, d = d)
  check_spd(ug = ug)
})

test_that("selective gram schmidt actually selects", {
  p <- 5
  d <- 0.25

  span <- matrix(ncol = p, nrow = p, stats::rnorm(p^2))

  madj <- upper.tri(matrix(NA, p, p))
  madj[madj] <- stats::rbinom(n = (p - 1) * p / 2, size = 1, prob = d)
  madj <- madj + t(madj)

  span_ort <- matrix(.C(
    "gram_schmidt_sel",
    double(p * p),
    as.logical(madj),
    as.double(t(span)),
    as.integer(p)
  )[[1]],
  ncol = p,
  byrow = TRUE
  )

  span_dot_prod <- tcrossprod(span_ort)
  madj_learned <- (zapsmall(span_dot_prod) != 0) * 1
  diag(madj_learned) <- 0
  expect_equal(length(which((madj_learned - madj) != 0)), 0)
})

test_that("the graph structure is preserved", {
  p <- 50
  d <- 0.2

  expect_equal_ug <- function(m, ug) {
    madj <- igraph::as_adjacency_matrix(ug, sparse = FALSE)
    madj_learned <- zapsmall(m) != 0
    diag(madj_learned) <- FALSE
    expect_equal(length(which((madj - madj_learned) != 0)), 0)
  }

  ug <- rgraph(p = p, d = d)

  sample <- port(ug = ug)
  expect_equal_ug(m = sample[, , 1], ug = ug)

  sample <- diagdom(ug = ug)
  expect_equal_ug(m = sample[, , 1], ug = ug)

  sample <- port_chol(ug = ug)
  expect_equal_ug(m = sample[, , 1], ug = ug)
})
