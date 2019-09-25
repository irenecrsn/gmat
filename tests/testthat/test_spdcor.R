context("Basic properties of the sampled matrices")

f_test <- c(
  chol_mh,
  chol_iid,
  port,
  diagdom,
  port_chol
)
f_test_dag <- c(chol_mh, chol_iid)

f_test_ug <- c(port, diagdom, port_chol)

test_that("the size of the sample is correct", {
  N <- 10

  for (f in f_test) {
    sample <- f(N = N)
    expect_equal(dim(sample)[3], N)
  }
})

test_that("matrix dimension is correct", {
  p <- 5

  for (f in f_test) {
    sample <- f(p = p)
    expect_equal(dim(sample)[1], dim(sample)[2])
    expect_equal(dim(sample)[1], p)
  }
})


test_that("matrices are symmetric positive definite", {
  p <- 5
  d <- 0.25

  check_spd <- function(sample) {
    expect_equal(sample[, , 1], t(sample[, , 1]))
    expect_gt(min(eigen(sample[, , 1])$values), 0)
  }

  for (f in f_test) {
    # no zeros
    sample <- f(p = p)
    check_spd(sample)

    # with a percentage of zeros
    sample <- f(p = p, d = d)
    check_spd(sample)
  }

  # with predefined zero patterns
  dag <- rgraph(p = p, d = d, dag = TRUE)
  for (f in f_test_dag) {
    sample <- f(dag = dag)
    check_spd(sample)
  }
  dag <- rgraph(p = p, d = d, dag = TRUE, ordered = FALSE)
  for (f in f_test_dag) {
    sample <- f(dag = dag)
    check_spd(sample)
  }
  ug <- rgraph(p = p, d = d)
  for (f in f_test_ug) {
    sample <- f(ug = ug)
    check_spd(sample)
  }
})

test_that("matrices are of correlation", {
  p <- 5
  d <- 0.25

  check_cor <- function(sample) {
    p <- dim(sample)[1]
    for (i in 1:(p - 1)) {
      expect_equal(sample[i, i, 1], 1)
      for (j in (i + 1):p) {
        expect_gt(sample[i, j, 1], -1)
        expect_lt(sample[i, j, 1], 1)
      }
    }
  }

  for (f in f_test) {
    # no zeros
    sample <- f(p = p)
    check_cor(sample)

    # with a percentage of zeros
    sample <- f(p = p, d = d)
    check_cor(sample)
  }

  # with predefined zero patterns
  dag <- rgraph(p = p, d = d, dag = TRUE)
  for (f in f_test_dag) {
    sample <- f(dag = dag)
    check_cor(sample)
  }
  dag <- rgraph(p = p, d = d, dag = TRUE, ordered = FALSE)
  for (f in f_test_dag) {
    sample <- f(dag = dag)
    check_cor(sample)
  }
  ug <- rgraph(p = p, d = d)
  for (f in f_test_ug) {
    sample <- f(ug = ug)
    check_cor(sample)
  }
})
