#' Sample a random Gaussian Markov network 
#'
#' Samples the concentration/covariance matrix corresponding to an 
#' undirected graph by two different methods: partial orthogonalization
#' or diagonal dominance.
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
#' @param k real number greater than `1`, the desired condition
#'number of the resulting matrix
#' @param zapzeros convert to zero extremely low entries
#'
#' @return  A three-dimensional array of length `p*p*N`
#'
#' @export
rgmn <- function(N = 1,
                 ug = NULL,
                 p = 10,
                 d = 0.25,
                 method = "port",
                 rentries = runif,
                 k = NULL,
								 zapzeros = FALSE) 
{
  
  if (is.null(ug)) {
    ug <- igraph::sample_gnp(n = p, p = d)
  }
  
  if (method == "port") {
    sam <- .rgmn_port(N = N, rentries = rentries, ug = ug, zapzeros = zapzeros)
  } else if (method == "diagdom") {
    sam <- .rgmn_diagdom(N = N, rentries = rentries, ug = ug)
  }
  
  if (!is.null(k) && k >= 1) {
    sam <- .kcontrol(sam = sam, k = k)	
  }
  
  return(sam)
}

.kcontrol <- function(sam, k = 1) {

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

.rgmn_diagdom <- function(N, rentries, ug) {
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

#' @useDynLib gmat, .registration=TRUE
.rgmn_port <- function(N, rentries, ug, zapzeros) {

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



