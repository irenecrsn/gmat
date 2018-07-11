context("port and domdiag methods")

# TODO change array so that we can do for (m in sample)
# (that is, first dimension is the number of samples)

test_that("the size of the sample is correct", {
	N <- 20; p <- 10; d <- 0.15;
	
	sample <- port(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)

	sample <- diagdom(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)
})

test_that("the concentration matrix is symmetric positive definite", {
	N <- 20; p <- 10; d <- 0.15;

	sample <- port(N = N, p = p, d = d) 
	for (n in 1:N) {
		expect_true(Matrix::isSymmetric(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}

	sample <- diagdom(N = N, p = p, d = d)
	for (n in 1:N) {
		expect_true(Matrix::isSymmetric(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
})

test_that("the graph structure is preserved", {
	N <- 20; p <- 10; d <- 0.15;

	expect_equal_ug <- function(m, ug) {
		madj <- as_adjacency_matrix(ug)
		madj_learned <- zapsmall(m) != 0
		diag(madj_learned) <- FALSE
		expect_equal(length(which((madj - madj_learned) != 0)), 0)
	}

	ug <- igraph::sample_gnp(n = p, p = d)

	sample <- port(N = N, p = p, ug = ug, d = d)
	for (n in 1:N) {
		expect_equal_ug(m = sample[, , n], ug = ug)
	}

	sample <- diagdom(N = N, p = p, ug = ug, d = d)
	for (n in 1:N) {
		expect_equal_ug(m = sample[, , n], ug = ug)
	}
})




