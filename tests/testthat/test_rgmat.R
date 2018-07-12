context("port and domdiag methods")

test_that("the condition number is controlled", {
	N <- 20; p <- 10; d <- 0.15; k <- 5;

	sample <- diagdom(N = N, p = p, d = d, k = k)
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}
})

test_that("selective gram schmidt actually selects", {
	p <- 10; d <- 0.15;

	span <- matrix(ncol = p, nrow = p, rnorm(p^2))
	
	madj <- upper.tri(matrix(NA, p, p))
	madj[madj] <- rbinom(n = (p - 1)*p / 2, size = 1, prob = d)
	madj <- madj + t(madj)

	span_ort <- matrix(.C("gram_schmidt_sel", 
							  double(p * p),
							  as.logical(madj),
							  as.double(t(span)),
							  as.integer(p))[[1]], 
						   ncol = p,
						   byrow = TRUE)
	
	span_dot_prod <- tcrossprod(span_ort)
	madj_learned <- (zapsmall(span_dot_prod) != 0) * 1
	diag(madj_learned) <- 0
	expect_equal(length(which((madj_learned - madj) != 0)), 0)
})

test_that("the size of the sample is correct", {
	N <- 20; p <- 10; d <- 0.15;
	
	sample <- port(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)

	sample <- diagdom(N = N, p = p, d = d)
	expect_equal(dim(sample)[3], N)
})

test_that("matrices are symmetric positive definite", {
	N <- 20; p <- 10; d <- 0.15;

	sample <- port(N = N, p = p, d = d) 
	for (n in 1:N) {
		expect_true(isSymmetric(sample[, , n]))
		# here we do not test positive definiteness since 
		# sometimes condition numbers are very high
	}

	sample <- diagdom(N = N, p = p, d = d)
	for (n in 1:N) {
		expect_true(isSymmetric(sample[, , n]))
		expect_gt(min(eigen(sample[, , n])$values), 0)
	}
})

test_that("the graph structure is preserved", {
	N <- 20; p <- 10; d <- 0.15;

	expect_equal_ug <- function(m, ug) {
		madj <- igraph::as_adjacency_matrix(ug, sparse = FALSE)
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




