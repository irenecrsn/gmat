context("Sampling of Gaussian Bayesian networks via Cholesky decomposition")

test_that("the size of the sample is correct", {
	N <- 20; p <- 10;

	sample <- rgbn_chol(N = N, p = p)
	expect_equal(dim(sample)[3], N)

	sample <- rgbn_iid(N = N, p = p)
	expect_equal(dim(sample)[3], N)

	sample <- rgbn_polar(N = N, p = p)
	expect_equal(dim(sample)[3], N)
})

test_that("matrix dimension is correct (full matrix)", {
	N <- 20; p <- 10;

	sample <- rgbn_chol(N = N, p = p)
	expect_equal(dim(sample)[1], dim(sample)[2])
	expect_equal(dim(sample)[1], p)

	sample <- rgbn_iid(N = N, p = p)
	expect_equal(dim(sample)[1], dim(sample)[2])
	expect_equal(dim(sample)[1], p)

	sample <- rgbn_polar(N = N, p = p)
	expect_equal(dim(sample)[1], dim(sample)[2])
	expect_equal(dim(sample)[1], p)
})

test_that("matrix dimension is correct (vectorized)", {
	N <- 20; p <- 10;
	p_vectorized <- p*(p - 1)/2

	sample <- rgbn_chol(N = N, p = p, return.minvector = TRUE)
	expect_equal(dim(sample)[2], p_vectorized)

	sample <- rgbn_iid(N = N, p = p, return.minvector = TRUE)
	expect_equal(dim(sample)[2], p_vectorized)

	sample <- rgbn_polar(N = N, p = p, return.minvector = TRUE)
	expect_equal(dim(sample)[2], p_vectorized)
})

test_that("matrices are symmetric positive definite", {
	N <- 20; p <- 10; d <- 0.15;

	sample <- rgbn_chol(N = N, p = p) 
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}

	sample <- rgbn_iid(N = N, p = p)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	
	sample <- rgbn_polar(N = N, p = p)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
})

test_that("the dag structure is preserved", {
	N <- 20; p <- 10; d <- 0.15;

	expect_equal_dag <- function(m, dag) {
		madj <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
		madj_learned <- zapsmall(m) != 0
		diag(madj_learned) <- FALSE
		expect_equal(length(which((madj - madj_learned) != 0)), 0)
	}

	dag <- rgraph(p = p, d = d, dag = TRUE)

	sample <- rgbn_chol(N = N, dag = dag, add_no_chordal = FALSE)
	for (n in 1:N) {
		L <- t(chol(anti_t(sample[, , n])))
		U <- t(anti_t(L))
		expect_equal_dag(m = U, dag = dag)
	}

	sample <- rgbn_iid(N = N, dag = dag)
	for (n in 1:N) {
		L <- t(chol(anti_t(sample[, , n])))
		U <- t(anti_t(L))
		expect_equal_dag(m = U, dag = dag)
	}
	
	sample <- rgbn_polar(N = N, dag = dag)
	for (n in 1:N) {
		L <- t(chol(anti_t(sample[, , n])))
		U <- t(anti_t(L))
		expect_equal_dag(m = U, dag = dag)
	}
})

