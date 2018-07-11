#' Sample random covariance/concentration graph matrices. 
#'
#' Samples the concentration/covariance matrix corresponding to an 
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
#' @param method one of the following: `port` (partial orthogonalization), 
#' `diagdom` (diagonal dominance)
#' @param rentries the random number generator for the
#'non-zero entries (defaults to `runif`)
#' @param zapzeros convert to zero extremely low entries
#'
#' @return  A three-dimensional array of length `p*p*N`
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
#' @export
diagdom <- function(N = 1, ug = NULL, p = 5, d = 0.25, rentries = runif) {
  
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
	return (sam)
}

#' Control the condition number of a sample. 
#'
#' Adjust the matrix sample provided so that all its elements have the desired
#' condition number
#'
#' @param sam matrix sample already generated
#' @param k real number greater than `1`, the desired condition
#'number of the resulting matrix
#'
#' @return  A three-dimensional array of length `p*p*N` with adjusted condition
#' number
#'
#' @export
kcontrol <- function(sam, k = 1) {

	N <- dim(sam)[3]
	p <- dim(sam)[1]

	for (n in 1:N) {
  		eig_val <- eigen(sam[, , n])$values
  		if (min(eig_val)<0){
			warning("matrix is not positive definite, return
					the original matrix")
  		} else {
  			delta <- (max(eig_val) - k*min(eig_val)) / (k - 1)
  			sam[, , n] <- sam[, , n] + 
				diag(x = delta, nrow = p)
		}
	}
  	return (sam)
}

