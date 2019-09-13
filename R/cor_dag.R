#' Simulation of correlation matrices
#'
#' Sample correlation matrices, possibly with a zero pattern in its Cholesky
#' decomposition constrained by an acyclic digraph.
#'
#' @name dag-constrained correlation matrices
#'
#' @rdname cor_dag
#'
#' @param N Number of samples.
#' @param p Matrix dimension. Ignored if `dag` is provided.
#' @param d Number in `[0,1]`, the proportion of non-zero
#' entries in the Cholesky factor of the sampled matrices.
#' Ignored if `dag` is provided.
#' @param dag An
#' igraph acyclic
#' digraph specifying the zero pattern in the upper Cholesky
#' factor of the sampled matrices. Nodes must be in ancestral
#' order, with the first one having no parents.
#' @param ... Additional parameters for [mh_u()].
#'
#' @details Function [chol_mh()] uses the method described in
#' Córdoba et al.  (2018) and implemented in [mh_u()], based
#' on a Metropolis-Hastings algorithm over the upper Cholesky
#' factorization.
#'
#' @return A three-dimensional array of length `p x p x N`.
#'
#' @references Córdoba I., Varando G., Bielza C., Larrañaga P. A fast
#' Metropolis-Hastings method for generating random correlation matrices. _Lecture Notes in
#' Computer Science_ (IDEAL 2018), vol 11314, pp. 117-124, 2018.
#'
#' @examples
#' ## Cholesky sampling via Metropolis-Hastings
#' # Generate a full matrix (default behaviour)
#' chol_mh()
#'
#' # Generate a matrix with a percentage of zeros
#' chol_mh(d = 0.5)
#'
#' # Generate a random acyclic digraph structure
#' dag <- rgraph(p = 3, d = 0.5, dag = TRUE)
#' igraph::print.igraph(dag)
#'
#' # Generate a matrix complying with the predefined zero pattern
#' chol_mh(dag = dag)
#' @export
chol_mh <- function(N = 1,
                    p = 3,
                    d = 1,
                    dag = NULL,
                    ...) {
  if (is.null(dag) == TRUE & d != 1) {
    # We generate the dag if a zero pattern is requested
    dag <- rgraph(p = p, d = d, dag = TRUE)
  }
  if (is.null(dag) == FALSE) {
    U <- mh_u(N = N, dag = dag, ...)
  } else {
    U <- mh_u(N = N, p = p, ...)
  }
  vC <- apply(U, MARGIN = 3, tcrossprod)
  C <- array(data = vC, dim = dim(U))

  return(C)
}

#' @rdname cor_dag
#'
#' @details The entries in the upper Cholesky factor are sampled i.i.d. by
#' function [chol_iid()], following Kalisch and Buhlmann (2007).
#'
#' @references Kalisch, M., Buhlmann, P. Estimating high-dimensional directed
#' acyclic graphs with the PC-algorithm, _Journal of Machine Learning Research_,
#' 8:613-636, 2007.
#'
#' @examples
#' ## Cholesky sampling via i.i.d. Cholesky factor
#' # Generate a full matrix (default behaviour)
#' chol_iid()
#'
#' # Generate a matrix with a percentage of zeros
#' chol_iid(d = 0.5)
#'
#' # Generate a matrix complying with the predefined zero pattern
#' igraph::print.igraph(dag)
#' chol_iid(dag = dag)
#' @export
chol_iid <- function(N = 1,
                     p = 3,
                     d = 1,
                     dag = NULL) {

  # We generate the dag if a zero pattern is requested
  if (is.null(dag) == TRUE & d != 1) {
    dag <- rgraph(p = p, d = d, dag = TRUE)
  }
  if (is.null(dag) == FALSE) {
    p <- length(igraph::V(dag))
    L_init <- t(igraph::as_adjacency_matrix(dag, sparse = FALSE))
    n_edges <- length(igraph::E(dag))
  } else {
    L_init <- matrix(nrow = p, ncol = p, data = 0)
    L_init[lower.tri(L_init)] <- 1
    n_edges <- p * (p - 1) / 2
  }
  R <- array(dim = c(p, p, N))

  for (n in 1:N) {
    L <- L_init
    L[L != 0] <- -runif(n = n_edges, min = 0.1, max = 1)
    diag(L) <- 1
    D <- diag(x = runif(p, 0.1, 1))
    Omega <- t(L) %*% solve(D) %*% L
    R[, , n] <- stats::cov2cor(Omega)
  }

  return(R)
}

#' Upper Cholesky factor sampling using Metropolis-Hastings
#'
#' Metropolis-Hasting algorithms to sample the upper Cholesky factor, using
#' positive hemispheres of different dimensions. A zero pattern may be specified
#' using an acyclic digraph.
#'
#' @name metropolis-hastings sampling
#'
#' @rdname mh
#'
#' @param N Number of samples.
#' @param p Dimension of the upper Cholesky factor.
#' @param dag An
#' igraph acyclic
#' digraph specifying the zero pattern in the upper Cholesky
#' factor of the sampled matrices. Nodes must be in ancestral
#' order, with the first one having no parents.
#' @param ... Additional parameters for [mh_sphere()].
#'
#' @author Gherardo Varando \email{gherardo.varando@math.ku.dk}
#'
#' @details Function [mh_u()] returns a sample of `N` upper Cholesky factors whose rows have
#' been generated using [mh_sphere()]. The dimensions of the hemispheres used to sample vary
#' depending both on the row number of the Cholesky factor, and whether there is a zero pattern
#' specified by `dag`.
#'
#' @examples
#' ## Upper Cholesky factor sampling
#' # Generate a random acyclic digraph
#' dag <- rgraph(p = 3, d = 0.5, dag = TRUE)
#' igraph::print.igraph(dag)
#'
#' # Generate an upper Cholesky factor complying with such zero pattern
#' mh_u(dag = dag)
#' # We may also generate it with no zero pattern (full upper triangular)
#' mh_u()
#' @export
mh_u <- function(N = 1,
                 p = 3,
                 dag = NULL,
                 ...) {
  if (is.null(dag) == FALSE) {
    p <- length(igraph::V(dag))
    dag_topo_sort <- as.numeric(igraph::topological.sort(dag))
    inv <- order(dag_topo_sort)
    u <- igraph::as_adjacency_matrix(dag, sparse = FALSE)[dag_topo_sort, dag_topo_sort]
    diag(u) <- 1
  }

  U <- array(dim = c(p, p, N), data = 0)
  U[p, p, 1:N] <- 1

  if (is.null(dag) == TRUE) {
    for (i in 1:(p - 1)) {
      su <- mh_sphere(N = N, k = p - i + 1, i = i, ...)
      U[i, i:p, 1:N] <- t(su)
    }
  } else {
    U[, , 1] <- u
    ch <- igraph::degree(dag, mode = "out")[dag_topo_sort]
    pa <- igraph::degree(dag, mode = "in")[dag_topo_sort]
    for (j in 1:(p - 1)) {
      su <- mh_sphere(N = N, k = ch[j] + 1, i = pa[j] + 1, ...)
      U[j, U[j, , 1] > 0, 1:N] <- t(su)
    }
    U <- array(data = U[inv, inv, ], dim = dim(U))
  }
  return(U)
}


#' @rdname mh
#'
#' @param k Dimension of the hemisphere from which the sample is taken.
#' @param i Integer, power of the first coordinate in the density.
#' @param h Heating phase size.
#' @param eps Perturbation variance.
#'
#' @details The details of the algorithm implemented by [mh_sphere()] can be found in the
#' paper Córdoba et al. (2018), including a discussion on
#' theoretical convergence and numerical experiments for
#' choosing its hyper parameters `h` and `eps`.
#'
#' @references Córdoba I., Varando G., Bielza C., Larrañaga P. A fast
#' Metropolis-Hastings method for generating random correlation matrices. _Lecture Notes in
#' Computer Science_ (IDEAL 2018), vol 11314, pp. 117-124, 2018.
#' @examples
#' ## Hemisphere sampling
#' # 3D hemisphere from a density proportional to the square of the first coordinate
#' mh_sphere(N = 4, k = 3, i = 2)
#' @export
mh_sphere <-
  function(N = 1,
             k,
             i = 1,
             h = 100,
             eps = 0.01) {
    Tot <- h + N # total number of iteration of MH
    Sample <- matrix(nrow = Tot, ncol = k) # obj initialization
    Sample[1, ] <- rnorm(n = k, mean = 0, sd = 1) # first point
    Sample[1, 1] <-
      abs(Sample[1, 1]) # absolute value first component (has to be positive)
    Sample[1, ] <-
      Sample[1, ] / sqrt(sum(Sample[1, ]^2)) # normalization
    for (j in 2:Tot) {
      prop <-
        Sample[j - 1, ] + rnorm(n = k, mean = 0, sd = eps) # perturbate previous sample
      prop <- prop / sqrt(sum(prop^2)) # normalize proposed
      if ((prop[1] > 0) &&
        (log(runif(1)) <= i * log((prop[1])) - i * log(Sample[j - 1, 1]))) {
        Sample[j, ] <- prop
      } else {
        Sample[j, ] <- Sample[j - 1, ]
      }
    }
    Sample <- Sample[(h + 1):Tot, ]

    return(Sample)
  }
