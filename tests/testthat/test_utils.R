context("Private utility functions")

test_that("projection is orthogonal", {
	p <- 10
	v <- rnorm(p)
	u <- rnorm(p)

	# Component along <u>
	v_proj_u <- .C("proj_ort", double(p), as.double(v),
				   as.double(u), as.integer(p))[[1]]

	v_orth_u <- v - v_proj_u # Component orthogonal to <u>
	expect_equal(sum(v_orth_u*u), 0)
	expect_equal(sum(v_proj_u*v_orth_u), 0)
})

test_that("gram schmidt process yields orthogonal base", {
	p <- 10
	span <- matrix(ncol = p, nrow = p, rnorm(p^2))
	
	# TODO don't use this trick and allow to call directly
	# "gram_schmidt"
	madj <- diag(1, nrow = p)

	span_ort <- matrix(.C("gram_schmidt_sel", 
							  double(p * p),
							  as.logical(madj),
							  as.double(t(span)),
							  as.integer(p))[[1]], 
						   ncol = p,
						   byrow = TRUE)
	span_dot_prod <- tcrossprod(span_ort)
	expect_true(isDiagonal(zapsmall(span_dot_prod)))
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

