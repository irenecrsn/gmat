context("Utility functions")

test_that("random dags are actually dags", {
  p <- 10
  d <- 0.5

  dag <- rgraph(p = p, d = d, dag = TRUE)
  expect_true(igraph::is_dag(dag))

  # Two random order checks
  dag <- rgraph(p = p, d = d, dag = TRUE, ordered = FALSE)
  expect_true(igraph::is_dag(dag))

  dag <- rgraph(p = p, d = d, dag = TRUE, ordered = FALSE)
  expect_true(igraph::is_dag(dag))
})

test_that("random dags follow the natural ancestral order by default", {
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

test_that("upper Cholesky factor is of Cholesky", {
  m <- port()[, , 1]
  u <- uchol(m)
  expect_equal(m, tcrossprod(u))
})

test_that("the dag orientation of an ug is actually a dag", {
  p <- 10
  d <- 0.5

  ug <- rgraph(p = p, d = d)
  dag <- ug_to_dag(ug = ug)
  expect_true(igraph::is_dag(dag))
})

test_that("the skeleton of the oriented dag is chordal", {
  p <- 10
  d <- 0.5

  ug <- rgraph(p = p, d = d)

  # We force a non chordal graph
  while (igraph::is_chordal(ug)$chordal == TRUE) {
    ug <- rgraph(p = p, d = d)
  }

  dag <- ug_to_dag(ug = ug)

  expect_true(igraph::is_chordal(igraph::as.undirected(dag))$chordal)
})

test_that("the skeleton of the oriented dag contains the original ug, keeping
the order", {
  p <- 10
  d <- 0.5

  ug <- rgraph(p = p, d = d)
  dag <- ug_to_dag(ug = ug)
  ug_cover <- igraph::as.undirected(dag)

  # This forces to keep the order
  domains <- as.list(1:p)

  expect_true(igraph::is_subgraph_isomorphic_to(ug, ug_cover, domains = domains))
})


test_that("the oriented dag does not contain v-structures", {
  p <- 10
  d <- 0.5

  ug <- rgraph(p = p, d = d)
  dag <- ug_to_dag(ug = ug)

  v_struct <- igraph::make_directed_graph(edges = c(1, 3, 2, 3))

  expect_false(igraph::is_subgraph_isomorphic_to(v_struct, dag, induced = TRUE))
})
