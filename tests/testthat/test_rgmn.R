context("Sampling GMNs: generalities")

# TODO change array so that we can do for (m in sample)
# (that is, first dimension is the number of samples)

test_that("the size of the sample is correct", {
	N <- 20; p <- 10; d <- 0.15;
	
	sample <- rgmn(N = N, p = p, d = d, method = "sqrt")
	expect_equal(dim(sample)[3], N)

	sample <- rgmn(N = N, p = p, d = d, method = "domdiag")
	expect_equal(dim(sample)[3], N)
})

test_that("the concentration matrix is symmetric positive definite", {
	N <- 20; p <- 10; d <- 0.15;

	sample <- rgmn(N = N, p = p, d = d, method = "sqrt", rentries = rnorm, zapzeros = TRUE)
	for (n in 1:N) {
		expect_true(isSymmetric(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}

	sample <- rgmn(N = N, p = p, d = d, method = "domdiag")
	for (n in 1:N) {
		expect_true(isSymmetric(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
})

test_that("the condition number is controlled", {
	N <- 20; p <- 10; d <- 0.15; k <- 5;

	sample <- rgmn(N = N, p = p, d = d, k = k, 
				   method = "sqrt")
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}

	sample <- rgmn(N = N, p = p, d = d, k = k, 
				   method = "domdiag")
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
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

	ug <- .rgraph(p = p, d = d)

	sample <- rgmn(N = N, p = p, d = d, method = "sqrt", 
				   ug = ug)
	for (n in 1:N) {
		expect_equal_ug(m = sample[, , n], ug = ug)
	}

	sample <- rgmn(N = N, p = p, d = d, method = "domdiag", 
				   ug = ug)
	for (n in 1:N) {
		expect_equal_ug(m = sample[, , n], ug = ug)
	}
})




