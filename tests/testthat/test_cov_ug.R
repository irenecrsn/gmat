context("SPD matrices, possibly with a zero pattern")

test_that("the size of the sample is correct", {
	N <- 10; p <- 5; d <- 0.25;

	check_sample_size <- function(N, ...) {
		sample <- port(N = N, ...)
		expect_equal(dim(sample)[3], N)
		sample <- diagdom(N = N, ...)
		expect_equal(dim(sample)[3], N)
	}

	# no zeros
	check_sample_size(N = N, p = p)

	# with a percentage of zeros
	check_sample_size(N = N, p = p, d = d)

	# with a predefined pattern of zeros
	ug <- rgraph(p = p, d = d)
	check_sample_size(N = N, ug = ug)
})

test_that("matrices are symmetric positive definite", {
	N <- 10; p <- 5; d <- 0.25;

	check_spd <- function(N = N, ...) {

		sample <- port(N = N, ...) 
		for (n in 1:N) {
			expect_equal(sample[, , n], t(sample[, , n]))
			# here we do not test positive definiteness since 
			# sometimes condition numbers are very high
		}
		sample <- diagdom(N = N, ...)
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
	ug <- rgraph(p = p, d = d)
	check_spd(N = N, ug = ug)
})

test_that("selective gram schmidt actually selects", {
	p <- 5; d <- 0.25;

	span <- matrix(ncol = p, nrow = p, stats::rnorm(p^2))
	
	madj <- upper.tri(matrix(NA, p, p))
	madj[madj] <- stats::rbinom(n = (p - 1)*p / 2, size = 1, prob = d)
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

test_that("the graph structure is preserved", {
	N <- 10; p <- 5; d <- 0.25;

	expect_equal_ug <- function(m, ug) {
		madj <- igraph::as_adjacency_matrix(ug, sparse = FALSE)
		madj_learned <- zapsmall(m) != 0
		diag(madj_learned) <- FALSE
		expect_equal(length(which((madj - madj_learned) != 0)), 0)
	}

	ug <- rgraph(p = p, d = d)

	sample <- port(N = N, ug = ug)
	for (n in 1:N) {
		expect_equal_ug(m = sample[, , n], ug = ug)
	}

	sample <- diagdom(N = N, ug = ug)
	for (n in 1:N) {
		expect_equal_ug(m = sample[, , n], ug = ug)
	}
})

test_that("the condition number is controlled", {
	N <- 10; p <- 5; d <- 0.25; k <- 5;

	# no zeros
	sample <- diagdom(N = N, p = p, k = k)
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}
	# with a percentage of zeros
	sample <- diagdom(N = N, p = p, d = d, k = k)
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}
	# with a predefined zero pattern 
	ug <- rgraph(p = p, d = d)
	sample <- diagdom(N = N, ug = ug, k = k)
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}
})




