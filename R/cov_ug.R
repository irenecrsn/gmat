#' Simulation of covariance matrices.
#'
#' Sample covariance matrices, possibly with a zero pattern constrained by an
#' undirected graph.
#'
#' @name ug-constrained covariance matrices
#' @rdname cov_ug
#'
#' @param N Number of samples.
#' @param p Matrix dimension. Ignored if `ug` is provided.
#' @param d Number in `[0,1]`, the proportion of non-zero
#' entries in the sampled matrices. Ignored if `ug` is provided.
#' @param ug An [igraph](https://CRAN.R-project.org/package=igraph) undirected graph specifying the zero pattern in the sampled matrices.
#' @param zapzeros Boolean, convert to zero extremely low entries? Defaults to `TRUE`.
#'
#' @details Function [port()] uses the method described in
#' Córdoba et al. (2018). In summary, it consists on generating a random
#' matrix `Q` and performing row-wise orthogonalization such that if `i` and `j`
#' are not adjacent in `ug`, then the rows corresponding to such indices are
#' orthogonalized, without violating previous orthogonalizations and without
#' introducing unwanted independences. The resulting matrix after the process
#' has finished is the cross product of `Q`.
#'
#' @return  A three-dimensional array of length `p x p x N`
#'
#' @references Córdoba, I., Varando, G., Bielza, C. and Larrañaga, P. A partial
#' orthogonalization method for simulation covariance and concentration graph
#' matrices. _Proceedings of Machine Learning Research_ (PGM 2018), vol. 72, pp.
#' 61 - 72, 2018.
#'
#' @examples
#' ## Partial orthogonalization
#' # Generate a full matrix (default behaviour)
#' port()
#' 
#' # Generate a matrix with a percentage of zeros
#' port(d = 0.5)
#' port(d = 0.5, zapzeros = FALSE) # no zero zap
#' 
#' # Generate a random undirected graph structure
#' ug <- rgraph(p = 3, d = 0.5)
#' igraph::print.igraph(ug)
#' 
#' # Generate a matrix complying with the predefined zero pattern
#' port(ug = ug)
#' port(ug = ug, zapzeros = FALSE) # no zero zap
#' @useDynLib gmat
#' @export
port <- function(N = 1, p = 3, d = 1, ug = NULL, zapzeros = TRUE) {
  if (is.null(ug) == TRUE & d != 1) {
    ug <- rgraph(p = p, d = d)
  }
  if (is.null(ug) == FALSE) {
    p <- length(igraph::V(ug))

    sam <- array(dim = c(p, p, N), data = 0)
    madj <- igraph::as_adjacency_matrix(ug,
      type = "both",
      sparse = FALSE
    )
    for (n in 1:N) {
      sam[, , n] <- matrix(nrow = p, ncol = p, data = runif(p^2))
      sam[, , n] <- matrix(.C(
        "gram_schmidt_sel",
        double(p * p),
        as.logical(madj),
        as.double(t(sam[, , n])),
        as.integer(p)
      )[[1]],
      ncol = p,
      byrow = TRUE
      )
      sam[, , n] <- tcrossprod(sam[, , n])

      if (zapzeros == TRUE) {
        sam[, , n] <- zapsmall(sam[, , n])
      }
    }
  } else {
    sam <- array(dim = c(p, p, N), data = runif(p * p * N))
    for (n in 1:N) {
      sam[, , n] <- tcrossprod(sam[, , n])
    }
  }

  return(sam)
}

#' @rdname cov_ug
#'
#' @details We also provide an implementation of the most commonly used in the
#' literature [diagdom()]. By contrast, this method produces a random matrix `M`
#' with zeros corresponding to missing edges in `ug`, and then enforces a
#' dominant diagonal to ensure positive definiteness. Matrices produced by
#' `diagdom` usually are better conditioned than those by `port`; however, they
#' typically suffer from small off-diagonal entries, which can compromise model
#' validation in Gaussian graphical models. This is avoided by `port`.
#'
#' @examples
#' ## Diagonal dominance
#' # Generate a full matrix (default behaviour)
#' diagdom()
#' 
#' # Generate a matrix with a percentage of zeros
#' diagdom(d = 0.5)
#' 
#' # Generate a matrix complying with the predefined zero pattern
#' igraph::print.igraph(ug)
#' diagdom(ug = ug)
#' @export
diagdom <- function(N = 1, p = 3, d = 1, ug = NULL) {

  # We generated the ug if a zero pattern is requested
  if (is.null(ug) == TRUE & d != 1) {
    ug <- rgraph(p = p, d = d)
  }
  if (is.null(ug) == FALSE) {
    p <- length(igraph::V(ug))
    edges <- igraph::as_edgelist(ug)

    sam <- array(dim = c(p, p, N), data = 0)
    ned <- nrow(edges)
    if (ned > 0) {
      for (i in 1:ned) {
        sam[edges[i, 1], edges[i, 2], ] <- sam[
          edges[i, 2],
          edges[i, 1],
        ] <-
          runif(N)
      }
    }
  } else {
    # Full covariance matrix
    sam <- array(dim = c(p, p, N))
    for (i in 1:p) {
      for (j in 1:p) {
        sam[i, j, ] <- sam[j, i, ] <- runif(N)
      }
    }
  }

  for (i in 1:p) {
    sam[i, i, ] <- abs(runif(N))
  }
  mdiag <- apply(sam,
    MARGIN = c(1, 3),
    FUN = function(row) {
      return(sum(abs(row)))
    }
  )
  sam <- sam + array(dim = dim(sam), data = apply(
    X = mdiag, MARGIN
    = 2, FUN = diag,
    nrow = p
  ))
  return(sam)
}
