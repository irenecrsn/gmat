context("Private utility functions")

test_that("random dags are actually dags", {
	p <- 10

	dag <- .rgraph(p, 0.5, dag = TRUE)
	expect_true(is_dag(dag))
})

test_that("random dags follow the natural ancestral order", {
	p <- 10
	
	dag <- .rgraph(p, 0.5, dag = TRUE)
	dag_am <- as_adjacency_matrix(dag)

	expect_true(isTriangular(dag_am, upper = TRUE))
})

test_that("random ugs are actually ugs", {
	p <- 10

	ug <- .rgraph(p, 0.5)
	ug_am <- as_adjacency_matrix(ug)

	expect_true(isSymmetric(ug_am))
	expect_identical(diag(ug_am), rep(0, p))
})



