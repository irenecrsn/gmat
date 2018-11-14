#' Sample a random Gaussian Bayesian network using different parametrizations
#' for the Cholesky decomposition. 
#' 
#' @name directed graph matrix simulation
#'
#' @rdname gbn-sim
#'
#' @param N Number of samples
#' @param p Number of variables
#' @param dag directed chordal acyclic graph, use the igraph package graph class 
#' @param return.minvector logical, if TRUE the minimimal vector representation
#' is returned (useful to plot in the elliptope)
#' @param add_no_chordal logical, if TRUE when the dag provided is not chordal,
#' a fill-in is computed, in order to ensure uniform distribution
#' @param ... additional parameters
#'
#' @details Function `rgbn_chol` uses the method described in Córdoba et al.
#' (2018), based on a Metropolis-Hastings algorithm over the upper Cholesky
#' factorization. 
#'
#' @return A three-dimensional array of length `p*p*N`
#'
#' @references Córdoba I., Varando G., Bielza C., Larrañaga P. A fast
#' Metropolis-Hastings method for generating random correlation matrices. Intelligent Data
#' Engineering and Automated Learning – IDEAL 2018. Lecture Notes in
#' Computer Science, vol 11314, pp. 117-124, 2018. 
#'
#' @export
rgbn_chol <- function(N = 1,
                  p = 10,
				  dag = NULL,
                  return.minvector = FALSE,
				  add_no_chordal = TRUE,
                  ...) {
	
	# Uniform sampling of chordal DAG
	if (is.null(dag) == FALSE & add_no_chordal == TRUE) {
   		isCh <- igraph::is_chordal(dag, fillin = TRUE)
   		
		if (isCh$chordal == FALSE){
    		dag <- igraph::add_edges(dag, edges = isCh$fillin)
   		}
	}
  sU <- mh_full(N = N, dag = dag, ...)
  vsC <- apply(sU, MARGIN = 3, function(U)
    return(U %*% t(U)))
  sC <- array(data = vsC, dim = dim(sU))
  if (return.minvector) {
    mv <- apply(sC, MARGIN = 3, function(m) {
      return(m[upper.tri(m)])
    })
    return(t(mv))
  } else{
    return(sC)
  }
}

#' @rdname gbn-sim
#'
#' @details The entries in the upper Cholesky factor are sampled i.i.d. by
#' function `rgbn_iid`. 
#'
#' @export
rgbn_iid <- function(N = 1,
				 p = 10,
				 dag = NULL,
				return.minvector = FALSE) 
{  	
	if (is.null(dag) == FALSE) {
		p <- length(igraph::V(dag))
		L_init <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
	} else {
		L_init <- matrix(nrow = p, ncol = p, data = 0)
		L_init[upper.tri(L_init)] <- 1
	}
	diag(L_init) <- 1
	R <- array(dim = c(p, p, N))

	for (n in 1:N) {
		L <- L_init
		D <- diag(x = runif(p, 0.1, 1))
		Omega <- t(L) %*% solve(D) %*% L 
		R[, , n] <- cov2cor(Omega) 
	}

  	if (return.minvector == TRUE) {
    	mv <- apply(R, MARGIN = 3, function(m) {
      		return(m[upper.tri(m)])
    	})
    	return(t(mv))
  	} else {
    	return(R)
  	}
}

#' @rdname gbn-sim
#'
#' @param comp String one of "numeric" or "recursive", indicating the
#' computational method to use for sampling the angles for "unifconc" method
#' 
#' @details Function `rgbn_polar` reparametrizes the Cholesky factor following
#' the approach by Pourahmadi and Wang (2015).
#' 
#' @references Pourahmadi, M., Wang, X. Distribution of random correlation matrices:
#' Hyperspherical parameterization of the Cholesky factor, Statistics &
#' Probability Letters, 106:5-12, 2015.
#'
#' @export
rgbn_polar <- function(N = 1,
				 p = 10,
                 comp = 'numeric',
				 dag = NULL,
				return.minvector = FALSE) 
{  	
	if (is.null(dag) == FALSE) {
		p <- length(igraph::V(dag))
		L_init <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
	} else {
		L_init <- matrix(nrow = p, ncol = p, data = 0)
		L_init[upper.tri(L_init)] <- 1
	}
	diag(L_init) <- 1
	R <- array(dim = c(p, p, N))

	for (n in 1:N) {
		U <- .rcoef_polar(p = p, method = comp, L = L_init)
		R[, , n] <- U %*% t(U)
	}

  	if (return.minvector == TRUE) {
    	mv <- apply(R, MARGIN = 3, function(m) {
      	return(m[upper.tri(m)])
    	})
    	return(t(mv))
  	} else{
    	return(R)
  	}
}


#' Sample coefficient matrix with polar parametrization
#'
#' Sample a polar parametrization of the Cholesky factor in the LL' decomposition of
#' the concentration matrix.
#'
#' @rdname polar
#'
#' @param p positive integer, the dimension of the square matrix
#' @param method String, method to sample from `sin(x)^k`: `numeric` or
#' `recursive`
#' @param edges Matrix defining the edges of the Bayesian network (from, to).
#'
#' @details The mcoef factor is sampled such that mcoefmcoef' is uniformly
#'distributed among the correlation matrices.
#'
#' @return `mcoef` a lower triangular matrix such that mcoefmcoef' has 1 on diagonal.
#'
.rcoef_polar <- function(p = 100,
                         method = 'numeric',
						 L) 
{	
	theta <- matrix(nrow = p, ncol = p, data = 0)
	theta[lower.tri(theta)] <- pi/2

	for (icol in 1:(p - 1)) {
		for (irow in (icol + 1):p) {
			if (L[irow, icol] != 0) {
				theta[irow, icol] <- .rsin(n = 1, k = p - icol, method = method)
				L[irow, icol] <- cos(theta[irow, icol])
			} 
		}
		if (icol >= 2) {
			for (j in 1:(icol - 1)) {
				L[icol:p, icol] <- L[icol:p, icol] * sin(theta[icol:p, j])
			}
		}
	}
	L[p, p] <- prod(sin(theta[p, 1:(p - 1)]))
	return(L)
}


#' @rdname polar
#'
#' @param k exponent for `sin^k`
#' @param x value where `sin^k` will be calculated
.sin_k <- function(x, k) {
	return (sin(x)^k)
}

#' @rdname polar
#'
#' @importFrom stats integrate
.sin_k_cum <- function(x, k, method = "numeric") {
	
	if (x <= 0) {
		return (0)
	} 
	
	if (x >= pi) {
		return (1)
	}

	if (method == "numeric") {
		const <- integrate(.sin_k, lower = 0, upper = pi, k = k)$value
		return(integrate(.sin_k, lower = 0, upper = x, k = k)$value / const)
	}
	const <- .sin_int(pi, k)
	return(.sin_int(x, k) / const)
}

#' @rdname polar
#' 
#' @details recursive computation 
.sin_int <- function(x, k = 1){
	if (length(x)>1){
		return(sapply(x,.sin_int,k))
	}
	if (x <= 0) {
		x<-0
	}
	if (x > pi) {
		x<-pi
	}
	if (k < 0) return(0)
	if (k == 0) {
		return(x)
	} else if (k == 1) {
		return(1-cos(x))
	} else {
		return(-(1/k)*cos(x)*.sin_k(x, k-1) + ((k-1)/k)*.sin_int(x, k-2)  )
	}
}

#' @rdname polar
#'
#' @param n Number of samples to generate
#'
#' @importFrom stats uniroot
.rsin <- function(n, k = 1, method = 'numeric') {
	.sin_k_inv_unif <- function(x, u) {
		return (.sin_k_cum(x, k, method) - u)
	}

	.sin_k_invsampl <- function(u) {
		return(uniroot(.sin_k_inv_unif, u = u, interval = c(0, pi), 
					   extendInt = "upX")$root)
	}

	return(sapply(runif(n), .sin_k_invsampl))
}


#' Full Metropolis-Hastings algorithm for correlation matrices/
#' Gaussian Bayesian networks
#'
#' @param N Number of samples
#' @param p Number of variables in the saturated model. This parameter is ignored if `dag` is provided
#' @param dag Acyclic digraph of the GBN to sample
#' @param h Burn-in phase
#' @param eps Perturbation variance
#'
#' @export
mh_full <- function(N = 1,
                    p = 10,
					dag = NULL,
                    h = 100,
                    eps = 0.1) {


	if (is.null(dag) == FALSE) {

  		p <- length(igraph::V(dag)) 
  		u <- igraph::as_adjacency_matrix(dag, sparse = FALSE)
  		diag(u)<-1
	}

  	U <- array(dim = c(p, p, N), data = 0)
  	U[p, p, 1:N] <- 1

	if (is.null(dag) == TRUE) {
  		for (i in 1:(p - 1)) {
    		su <- mh_row(N = N,	p = p - i + 1, i = i, h = h, eps = eps)
    		U[i, i:p, 1:N] <- t(su)
  		}
	} else {
		U[, , 1] <- u
  		ch <- igraph::degree(dag, mode = "out")
  		pa <- igraph::degree(dag, mode = "in")
  		for (j in 1:(p-1)) {
    		su <- mh_row(N = N, p = ch[j] + 1, i = pa[j] + 1, h = h, eps=eps)
    		U[j, U[j, , 1] > 0, 1:N] <- t(su)
  		}
	}
  return(U)
}


#' sampling on sphere proportionally to a power of the first coordinate
#' 
#' @param N sample size
#' @param p dimension
#' @param i exponent of the density
#' @param eps Perturbation variance
#' @param returnAll Include in the output samples from the heat-in phase
#' @param h heating phase size
#' 
#'  Metropolis-Hasting algorithm to sample in the n dimensional semi-sphere (x_1>0)  
#' @export
mh_row <-
  function(N = 1,
           p,
           i = 1,
           h = 100,
           eps = 0.01,
           returnAll = FALSE) {
    Tot <- h + N #total number of iteration of MH
    Sample <- matrix(nrow = Tot, ncol = p) #obj initialization
    Sample[1, ] <- rnorm(n = p, mean = 0, sd = 1) #first point
    Sample[1, 1] <-
      abs(Sample[1, 1]) #absolute value first component (has to be positive)
    Sample[1, ] <-
      Sample[1, ] / sqrt(sum(Sample[1,] ^ 2)) #normalization
    for (j in 2:Tot) {
      prop <-
        Sample[j - 1, ] + rnorm(n = p, mean = 0, sd = eps) # perturbate previous sample
      prop <- prop / sqrt(sum(prop ^ 2)) #normalize proposed
      if ((prop[1] > 0) &&
          (log(runif(1)) <= i * log((prop[1])) - i * log(Sample[j - 1 , 1]))) {
        Sample[j, ] <- prop
      } else{
        Sample[j, ] <- Sample[j - 1 ,]
      }
    }
    if (returnAll == FALSE) {
      Sample <- Sample[(h + 1):Tot, ]
    }
    
    return(Sample)
  }
