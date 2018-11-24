context("Correlation matrices, possibly with a zero pattern on the Cholesky factor")

test_that("the size of the sample is correct", {
	N <- 10; p <- 5; d <- 0.25;

	check_sample_size <- function(N, ...) {
		sample <- chol_mh(N = N, ...)
		expect_equal(dim(sample)[3], N)
		sample <- chol_iid(N = N, ...)
		expect_equal(dim(sample)[3], N)
		sample <- chol_polar(N = N, ...)
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


test_that("matrices are symmetric positive definite", {
	N <- 10; p <- 5; d <- 0.25;

	check_spd <- function(N = N, ...) {
		sample <- chol_mh(N = N, ...) 
		for (n in 1:N) {
			expect_equal(sample[, , n], t(sample[, , n]))
			expect_gt(min(eigen(sample[, , n])$values), 0)
		}
		sample <- chol_iid(N = N, ...)
		for (n in 1:N) {
			expect_equal(sample[, , n], t(sample[, , n]))
			expect_gt(min(eigen(sample[, , n])$values), 0)
		}
		sample <- chol_polar(N = N, ...)
		for (n in 1:N) {
			expect_equal(sample[, , n], t(sample[, , n]))
			expect_gt(min(eigen(sample[, , n])$values), 0)
		}
	}

	# no zeros
	check_spd(N = N, p = p)

	# with a percentage of zeros
	check_spd(N = N, p = p, d = d)

	# with a predefined zero pattern
	dag <- rgraph(p = p, d = d, dag = TRUE)
	check_spd(N = N, dag = dag)
})

test_that("the dag structure is preserved", {
	N <- 10; p <- 5; d <- 0.25;

	expect_equal_dag <- function(m, dag) {
		L <- t(chol(anti_t(m)))
		U <- t(anti_t(L))
		madj <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
		madj_learned <- zapsmall(U) != 0
		diag(madj_learned) <- FALSE
		expect_equal(length(which((madj - madj_learned) != 0)), 0)
	}

	dag <- rgraph(p = p, d = d, dag = TRUE)

	sample <- chol_mh(N = N, dag = dag)
	for (n in 1:N) {
		expect_equal_dag(m = sample[, , n], dag = dag)
	}

	sample <- chol_iid(N = N, dag = dag)
	for (n in 1:N) {
		expect_equal_dag(m = sample[, , n], dag = dag)
	}
	
	sample <- chol_polar(N = N, dag = dag)
	for (n in 1:N) {
		expect_equal_dag(m = sample[, , n], dag = dag)
	}
})

