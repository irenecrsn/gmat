#' Sample random covariance/concentration graph matrices. 
#'
#' Samples the covariance/concentration matrix corresponding to an 
#' undirected graph by either partial orthogonalization or diagonal dominance.
#'
#' @name rgmat
#' @rdname rgmat
#'
#' @param N sample size
#' @param ug the undirected graph 
#' @param p Matrix dimension (number of nodes in `ug` if provided) 
#' @param d number in `[0,1]`, the proportion of non-zero
#'entries (if `ug` is not provided)
#' @param rentries the random number generator for the
#'non-zero entries (defaults to `runif`)
#' @param zapzeros convert to zero extremely low entries
#' 
#' @details Function `port` uses the method described in 
#' [https://arxiv.org/abs/1807.03090](https://arxiv.org/abs/1807.03090). In summary, it consists on generating a random
#' matrix `Q` and performing row-wise orthogonalization such that if `i` and `j`
#' are not adjacent in `ug`, then the rows corresponding to such indices are
#' orthogonalized, without violating previous orthogonalizations and without
#' introducing unwanted independences. The resulting matrix after the process
#' has finished is the cross product of `Q`.
#'
#' @return  A three-dimensional array of length `p*p*N`
#' 
#' @examples
#'
#' # Generate a random undirected graph structure
#' ug <- igraph::sample_gnp(n = 3, p = 0.25)
#'
#' # Generate 10 matrices complying with such random structure via 
#' # partial orthogonalization
#' gmat::port(N = 10, ug = ug)
#'
#' @useDynLib gmat, .registration=TRUE
#' @export
port <- function(N = 1, ug = NULL, p = 5, d = 0.25, rentries = runif, zapzeros = TRUE) {
  
  if (is.null(ug)) {
    ug <- igraph::sample_gnp(n = p, p = d)
  }
	p <- length(V(ug))

	sam <- array(dim = c(p, p, N), data = 0)

	madj <- igraph::as_adjacency_matrix(ug, type = "both", 
										sparse = FALSE)
	for (n in 1:N) {
	  sam[, , n] <- matrix(nrow = p, ncol = p, data = rentries(p^2))
	  sam[, , n] <- matrix(.C("gram_schmidt_sel", 
							  double(p * p),
							  as.logical(madj),
							  as.double(t(sam[, , n])),
							  as.integer(p))[[1]], 
						   ncol = p,
						   byrow = TRUE)
	  sam[, , n] <- tcrossprod(sam[, , n])
	  
	  if (zapzeros == TRUE) {
	  	sam[, , n] <- zapsmall(sam[, , n])
	  }
	}
	 
	return(sam)
}

#' @rdname rgmat
#'
#' @param k real number greater than `1`, the desired condition
#'	number of the matrices in the resulting sample 
#'
#' @details We also provide an implementation of the most commonly used in the
#' literature `diagdom`. By contrast, this method produces a random matrix `M`
#' with zeros corresponding to missing edges in `ug`, and then enforces a
#' dominant diagonal to ensure positive definiteness. Matrices produced by
#' `diagdom` usually are better conditioned than those by `port`; however, they
#' typically suffer from small off-diagonal entries, which can compromise model
#' validation in Gaussian graphical models. This is avoided by `port`.
#' 
#' @examples
#' # Generate 10 matrices complying with such random structure via 
#' # diagonal dominance 
#' gmat::diagdom(N = 10, ug = ug)
#'
#' @export
diagdom <- function(N = 1, ug = NULL, p = 5, d = 0.25, rentries = runif, k = NULL) {
  
  if (is.null(ug)) {
    ug <- igraph::sample_gnp(n = p, p = d)
  }
  p <- length(V(ug))
  edges <- as_edgelist(ug)

  sam <- array(dim = c(p, p, N), data = 0)
  ned <- nrow(edges)
  if (ned > 0) {
    for (i in 1:ned) {
      sam[edges[i, 1], edges[i, 2], ] <- sam[edges[i ,2],
										    edges[i, 1], ] <-
											    rentries(N) 
    }
  }
  for (i in 1:p) {
    sam[i, i, ] <- abs(rentries(N)) 
  }
  mdiag <- apply(sam, MARGIN = c(1, 3), 
                 FUN = function(row) {return(sum(abs(row)))})
  sam <- sam + array(dim = dim(sam), data = apply(X = mdiag, MARGIN
												  = 2, FUN = diag,
												  nrow = p))
	if (!is.null(k)) {
		for (n in 1:N) {
  			eig_val <- eigen(sam[, , n])$values
  			delta <- (max(eig_val) - k*min(eig_val)) / (k - 1)
  			sam[, , n] <- sam[, , n] + diag(x = delta, nrow = p)
		}
	}

	return (sam)
}

