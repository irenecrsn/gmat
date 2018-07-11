context("Condition number control")

test_that("the condition number is controlled", {
	N <- 20; p <- 10; d <- 0.15; k <- 5;

	sample <- port(N = N, p = p, d = d)
	sample <- kcontrol(sam = sample, k = k)
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}

	sample <- diagdom(N = N, p = p, d = d)
	sample <- kcontrol(sam = sample, k = k)
	for (n in 1:N) {
		expect_equal(kappa(sample[, , n], exact = TRUE), k)
	}
})

