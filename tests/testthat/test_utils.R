context("Utility functions")

test_that("random dags are actually dags", {
  p <- 10

  dag <- rgraph(p, 0.5, dag = TRUE)
  expect_true(igraph::is_dag(dag))
})

test_that("random dags follow the natural ancestral order", {
  p <- 10

  dag <- rgraph(p, 0.5, dag = TRUE)
  dag_am <- igraph::as_adjacency_matrix(dag, sparse = FALSE)

  # We simply check if the matrix is upper triangular
  expect_equal(sum(diag(dag_am)), 0)
  expect_equal(sum(dag_am[lower.tri(dag_am)]), 0)
})

test_that("random ugs are actually ugs", {
  p <- 10

  ug <- rgraph(p, 0.5)
  ug_am <- igraph::as_adjacency_matrix(ug, sparse = FALSE)

  expect_equal(ug_am, t(ug_am))
  expect_equal(sum(diag(ug_am)), 0)
})

test_that("sample vectorization works correctly", {
  N <- 10
  p <- 5
  p_vectorized <- p * (p - 1) / 2

  sample <- diagdom(N = N, p = p)
  sample <- vectorize(sample)
  expect_equal(length(dim(sample)), 2)
  expect_equal(dim(sample)[2], p_vectorized)
})

test_that("the condition number is correctly set", {
  N <- 10
  p <- 5
  k <- 5

  sample <- diagdom(N = N, p = p)
  sample <- set_cond_number(sample = sample, k = k)
  for (n in 1:N) {
    expect_equal(kappa(sample[, , n], exact = TRUE), k)
  }
})
