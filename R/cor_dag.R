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
#' Ignored if `dag` is provided. Ignored by [chol_polar()]
#' @param dag An
#' [igraph](https://CRAN.R-project.org/package=igraph) acyclic
#' digraph specifying the zero pattern in the upper Cholesky
#' factor of the sampled matrices. Nodes must be in ancestral
#' order, with the first one having no parents.
#' @param add_no_chordal Logical, if TRUE when the dag provided is not chordal,
#' a fill-in is computed, in order to ensure uniform distribution. Ignored if
#' `dag` or `d` are not provided. Defaults to FALSE.
#' @param ... additional parameters
#'
#' @details Function [chol_mh()] uses the method described in Córdoba et al.
#' (2018), based on a Metropolis-Hastings algorithm over the upper Cholesky
#' factorization.
#'
#' @return A three-dimensional array of length `p x p x N`
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
                    add_no_chordal = FALSE,
                    ...) {
  if (is.null(dag) == TRUE) {
    # We generate the dag if a zero pattern is requested
    if (d != 1) {
      dag <- rgraph(p = p, d = d, dag = TRUE)
    } else {
      sU <- .rcoef_mh(N = N, p = p, ...)
    }
  }
  if (is.null(dag) == FALSE) {
    # Uniform sampling of chordal DAG
    if (add_no_chordal == TRUE) {
      isCh <- igraph::is_chordal(dag, fillin = TRUE)

      if (isCh$chordal == FALSE) {
        dag <- igraph::add_edges(dag, edges = isCh$fillin)
      }
    }
    sU <- .rcoef_mh(N = N, dag = dag, ...)
  }
  vsC <- apply(sU, MARGIN = 3, tcrossprod)
  sC <- array(data = vsC, dim = dim(sU))

  return(sC)
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

#' @rdname cor_dag
#'
#' @param comp String one of "numeric" or "recursive", indicating the
#' computational method to use for sampling the angles for "unifconc" method
#'
#' @details Function [chol_polar()] reparametrizes the Cholesky factor following
#' the approach by Pourahmadi and Wang (2015).
#'
#' @references Pourahmadi, M., Wang, X. Distribution of random correlation matrices:
#' Hyperspherical parameterization of the Cholesky factor, _Statistics &
#' Probability Letters_, 106:5-12, 2015.
#'
#' @examples
#' ## Cholesky sampling via polar parametrization of the lower Cholesky factor
#' # Generate a full matrix (default behaviour)
#' chol_polar()
#' 
#' # Generate a matrix with a percentage of zeros
#' chol_polar(d = 0.5)
#' 
#' # Generate a matrix complying with the predefined zero pattern
#' igraph::print.igraph(dag)
#' chol_polar(dag = dag)
#' 
#' # Performance comparison of numeric vs recursive integral (full matrix)
#' system.time(chol_polar(N = 10, p = 5))
#' system.time(chol_polar(N = 10, p = 5, comp = "recursive"))
#' @export
chol_polar <- function(N = 1, p = 3,
                       d = 1,
                       dag = NULL,
                       comp = "numeric") {
  # We generate the dag if a zero pattern is requested
  if (is.null(dag) == TRUE & d != 1) {
    dag <- rgraph(p = p, d = d, dag = TRUE)
  }
  if (is.null(dag) == FALSE) {
    p <- length(igraph::V(dag))
    L_init <- t(igraph::as_adjacency_matrix(dag, sparse = FALSE))
  } else {
    L_init <- matrix(nrow = p, ncol = p, data = 0)
    L_init[lower.tri(L_init)] <- 1
  }
  diag(L_init) <- 1
  R <- array(dim = c(p, p, N))

  for (n in 1:N) {
    L <- .rcoef_polar(p = p, method = comp, L = L_init)
    R[, , n] <- crossprod(L)
  }

  return(R)
}

.rcoef_polar <- function(p, method, L) {
  theta <- matrix(nrow = p, ncol = p, data = 0)
  theta[lower.tri(theta)] <- pi / 2

  for (j in 1:(p - 1)) {
    for (i in (j + 1):p) {
      if (L[i, j] != 0) {
        theta[i, j] <- .rsin(n = 1, k = p - j, method = method)
        L[i, j] <- cos(theta[i, j])
      }
    }
    if (j >= 2) {
      for (k in 1:(j - 1)) {
        L[j:p, j] <- L[j:p, j] * sin(theta[j:p, k])
      }
    }
  }
  L[p, p] <- prod(sin(theta[p, 1:(p - 1)]))
  return(L)
}


.sin_k <- function(x, k) {
  return(sin(x)^k)
}

.sin_k_cum <- function(x, k, method = "numeric") {
  if (x <= 0) {
    return(0)
  }

  if (x >= pi) {
    return(1)
  }

  if (method == "numeric") {
    const <- stats::integrate(.sin_k, lower = 0, upper = pi, k = k)$value
    return(stats::integrate(.sin_k, lower = 0, upper = x, k = k)$value / const)
  } else {
    const <- .sin_int(pi, k)
    return(.sin_int(x, k) / const)
  }
}

.sin_int <- function(x, k = 1) {
  if (length(x) > 1) {
    return(sapply(x, .sin_int, k))
  }
  if (x <= 0) {
    x <- 0
  }
  if (x > pi) {
    x <- pi
  }
  if (k < 0) return(0)
  if (k == 0) {
    return(x)
  } else if (k == 1) {
    return(1 - cos(x))
  } else {
    return(-(1 / k) * cos(x) * .sin_k(x, k - 1) + ((k - 1) / k) * .sin_int(x, k - 2))
  }
}

.rsin <- function(n, k = 1, method = "numeric") {
  .sin_k_inv_unif <- function(x, u) {
    return(.sin_k_cum(x, k, method) - u)
  }

  .sin_k_invsampl <- function(u) {
    return(stats::uniroot(.sin_k_inv_unif,
      u = u, interval = c(0, pi),
      extendInt = "upX"
    )$root)
  }

  return(sapply(runif(n), .sin_k_invsampl))
}


.rcoef_mh <- function(N = 1,
                      p = 3,
                      dag = NULL,
                      h = 100,
                      eps = 0.1) {
  if (is.null(dag) == FALSE) {
    p <- length(igraph::V(dag))
    u <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
    diag(u) <- 1
  }

  U <- array(dim = c(p, p, N), data = 0)
  U[p, p, 1:N] <- 1

  if (is.null(dag) == TRUE) {
    for (i in 1:(p - 1)) {
      su <- mh_sphere(N = N, p = p - i + 1, i = i, h = h, eps = eps)
      U[i, i:p, 1:N] <- t(su)
    }
  } else {
    U[, , 1] <- u
    ch <- igraph::degree(dag, mode = "out")
    pa <- igraph::degree(dag, mode = "in")
    for (j in 1:(p - 1)) {
      su <- mh_sphere(N = N, p = ch[j] + 1, i = pa[j] + 1, h = h, eps = eps)
      U[j, U[j, , 1] > 0, 1:N] <- t(su)
    }
  }
  return(U)
}


#' Sampling on sphere proportionally to a power of the first coordinate
#'
#' Metropolis-Hasting algorithm to sample in the
#' `p`-dimensional hemisphere.
#'
#' @param N Number of samples.
#' @param p Dimension of the hemisphere.
#' @param i Integer, power of the first coordinate in the density.
#' @param eps Perturbation variance.
#' @param h Heating phase size
#'
#' @author Gherardo Varando \email{gherardo.varando@math.ku.dk}
#'
#' @details The details of the algorithm can be found in the paper Córdoba et al. (2018), including a discussion on theoretical convergence and numerical experiments for choosing its hyper parameters.
#'
#' @references Córdoba I., Varando G., Bielza C., Larrañaga P. A fast
#' Metropolis-Hastings method for generating random correlation matrices. _Lecture Notes in
#' Computer Science_ (IDEAL 2018), vol 11314, pp. 117-124, 2018.
#' @examples
#' mh_sphere(N = 4, p = 3, i = 2)
#' @export
mh_sphere <-
  function(N = 1,
             p,
             i = 1,
             h = 100,
             eps = 0.01) {
    Tot <- h + N # total number of iteration of MH
    Sample <- matrix(nrow = Tot, ncol = p) # obj initialization
    Sample[1, ] <- rnorm(n = p, mean = 0, sd = 1) # first point
    Sample[1, 1] <-
      abs(Sample[1, 1]) # absolute value first component (has to be positive)
    Sample[1, ] <-
      Sample[1, ] / sqrt(sum(Sample[1, ]^2)) # normalization
    for (j in 2:Tot) {
      prop <-
        Sample[j - 1, ] + rnorm(n = p, mean = 0, sd = eps) # perturbate previous sample
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
