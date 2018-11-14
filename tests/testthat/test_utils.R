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
	expect_identical(sum(diag(dag_am)), 0)
	expect_identical(sum(dag_am[lower.tri(dag_am)]), 0)
})

test_that("random ugs are actually ugs", {
	p <- 10

	ug <- rgraph(p, 0.5)
	ug_am <- igraph::as_adjacency_matrix(ug, sparse = FALSE)

	expect_equal(ug_am, t(ug_am))
	expect_identical(sum(diag(ug_am)), 0)
})



