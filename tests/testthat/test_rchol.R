context("Utility functions")

test_that("random dags are actually dags", {
	p <- 10

	dag <- rgraph(p, 0.5, dag = TRUE)
	expect_true(is_dag(dag))
})

test_that("random dags follow the natural ancestral order", {
	p <- 10
	
	dag <- rgraph(p, 0.5, dag = TRUE)
	dag_am <- igraph::as_adjacency_matrix(dag)

	expect_true(Matrix::isTriangular(dag_am, upper = TRUE))
})

test_that("random ugs are actually ugs", {
	p <- 10

	ug <- rgraph(p, 0.5)
	ug_am <- igraph::as_adjacency_matrix(ug)

	expect_true(Matrix::isSymmetric(ug_am))
	expect_identical(sum(Matrix::diag(ug_am)), 0)
})



