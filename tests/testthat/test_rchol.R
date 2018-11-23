context("SPD matrices, possibly with a zero pattern on the Cholesky factor")

test_that("the size of the sample is correct", {
	N <- 10; p <- 5; d <- 0.25;

	# no zeros
	sample <- chol_mh(N = N, p = p)
	expect_equal(dim(sample)[3], N)
	sample <- chol_iid(N = N, p = p)
	expect_equal(dim(sample)[3], N)
	sample <- chol_polar(N = N, p = p)
	expect_equal(dim(sample)[3], N)

	# with a percentage of zeros
	sample <- chol_mh(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)
	sample <- chol_iid(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)
	sample <- chol_polar(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)

	# with a predefined pattern of zeros
	dag <- rgraph(p = p, d = d, dag = TRUE)
	sample <- chol_mh(N = N, dag = dag)
	expect_equal(dim(sample)[3], N)
	sample <- chol_iid(N = N, dag = dag)
	expect_equal(dim(sample)[3], N)
	sample <- chol_polar(N = N, dag = dag)
	expect_equal(dim(sample)[3], N)
})

test_that("matrix dimension is correct (full matrix)", {
	N <- 10; p <- 5;

	sample <- chol_mh(N = N, p = p)
	expect_equal(dim(sample)[1], dim(sample)[2])
	expect_equal(dim(sample)[1], p)
	sample <- chol_iid(N = N, p = p)
	expect_equal(dim(sample)[1], dim(sample)[2])
	expect_equal(dim(sample)[1], p)
	sample <- chol_polar(N = N, p = p)
	expect_equal(dim(sample)[1], dim(sample)[2])
	expect_equal(dim(sample)[1], p)
})

test_that("matrix dimension is correct (vectorized)", {
	N <- 10; p <- 5;
	p_vectorized <- p*(p - 1)/2

	sample <- chol_mh(N = N, p = p, return.minvector = TRUE)
	expect_equal(dim(sample)[2], p_vectorized)

	sample <- chol_iid(N = N, p = p, return.minvector = TRUE)
	expect_equal(dim(sample)[2], p_vectorized)

	sample <- chol_polar(N = N, p = p, return.minvector = TRUE)
	expect_equal(dim(sample)[2], p_vectorized)
})

test_that("matrices are symmetric positive definite", {
	N <- 10; p <- 5; d <- 0.25;

	# no zeros
	sample <- chol_mh(N = N, p = p) 
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	sample <- chol_iid(N = N, p = p)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	sample <- chol_polar(N = N, p = p)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}

	# with a percentage of zeros
	sample <- chol_mh(N = N, p = p, d = d) 
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	sample <- chol_iid(N = N, p = p, d = d)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	sample <- chol_polar(N = N, p = p, d = d)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	
	# with a predefined zero pattern
	dag <- rgraph(p = p, d = d, dag = TRUE)
	sample <- chol_mh(N = N, dag = dag) 
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	sample <- chol_iid(N = N, dag = dag)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	sample <- chol_polar(N = N, dag = dag)
	for (n in 1:N) {
		expect_equal(sample[, , n], t(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
	
})

test_that("the dag structure is preserved", {
	N <- 10; p <- 5; d <- 0.25;

	expect_equal_dag <- function(m, dag) {
		madj <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
		madj_learned <- zapsmall(m) != 0
		diag(madj_learned) <- FALSE
		expect_equal(length(which((madj - madj_learned) != 0)), 0)
	}

	dag <- rgraph(p = p, d = d, dag = TRUE)

	sample <- chol_mh(N = N, dag = dag)
	for (n in 1:N) {
		L <- t(chol(anti_t(sample[, , n])))
		U <- t(anti_t(L))
		expect_equal_dag(m = U, dag = dag)
	}

	sample <- chol_iid(N = N, dag = dag)
	for (n in 1:N) {
		L <- t(chol(anti_t(sample[, , n])))
		U <- t(anti_t(L))
		expect_equal_dag(m = U, dag = dag)
	}
	
	sample <- chol_polar(N = N, dag = dag)
	for (n in 1:N) {
		L <- t(chol(anti_t(sample[, , n])))
		U <- t(anti_t(L))
		expect_equal_dag(m = U, dag = dag)
	}
})

