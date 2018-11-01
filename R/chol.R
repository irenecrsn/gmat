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
      prop <- prop / sqrt(sum(prop ^ 2)) #normalize propozed
      if ((prop[1] > 0) &&
          (log(runif(1)) <= i * log((prop[1])) - i * log(Sample[j - 1 , 1]))) {
        Sample[j, ] <- prop
      } else{
        Sample[j, ] <- Sample[j - 1 ,]
      }
    }
    if (!returnAll) {
      Sample <- Sample[(h + 1):Tot,]
    }
    
    return(Sample)
  }


mh_full <- function(N = 1,
                    p = 10,
                    h = 100,
                    eps = 0.1) {
  U <- array(dim = c(p, p, N), data = 0)
  U[p, p, 1:N] <- 1
  for (i in 1:(p - 1)) {
    su <- mh_row(
      N = N,
      p = p - i + 1 ,
      i = i ,
      h = h,
      eps = eps
    )
    U[i, i:p, 1:N] <- t(su)
  }
  return(U)
}

#' uchol
#'
#'
#' @param  N
#' @param return.minvector logical, if TRUE the minimimal vector representation
#' is returned (useful to plot in the elliptope)
#' @param ... additional parameters
#'
#' @export
uchol <- function(N = 1,
                  p = 10,
                  return.minvector = FALSE,
                  ...) {
  sU <- mh_full(N = N, p = p, ...)
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

#' Sample a random Gaussian Bayesian network
#' 
#' Samples a random Gaussian Bayesian network with different methods.
#'
#' @param p Positive integer, Number of dimension.
#' @param method String one of "iidcoef", "unifconc"
#' @param d number in `[0,1]`, the proportion of non-zero entries.
#' @param rmcoef function, the random number generator for the non-zero
#' entries of the coefficient matrix for the "iidcoef" method.
#' @param comp String one of "numeric" or "recursive", indicating the
#' computational method to use for sampling the angles for "unifconc" method
#' @param dag A graph object of class `igraph`.
#'
#' @return  A list with the following components:
#'   - param A list with the parameters for the concentration matrix. For 
#'   "iidcoef", these are the coefficient matrix and the conditional variance
#'   matrix. For "unifconc", the parameters are the angle matrix and the
#'   coefficient matrix.
#'   - mconc The sampled concentration matrix 
#'   - dag The resulting dag structure
#'   
#' @export
rgbn <- function(p = 10,
                 method = "iidcoef",
                 d = 0.25,
                 rmcoef = rnorm, 
                 comp = 'numeric',
				 dag = NULL) 
{  	
	if (is.null(dag)) {
		dag <- .rgraph(p, d, dag = TRUE)
	} 	
	edges <- as_edgelist(dag)

	if (method == "iidcoef") {
		param <- list()
		mcoef <- .rcoef_iid(p = p, rmcoef = rmcoef, edges = edges)
		mccov <- Diagonal(n = p, runif(p, 0.1, 1))
		mconc <- t(mcoef) %*% solve(mccov) %*% (mcoef) 
		param <- list(mcoef = mcoef, mccov = mccov)
	} else if (method == "unifconc") {
		param <- .rcoef_polar(p = p, method = comp, edges = edges)
		mconc <- param$mcoef %*% t(param$mcoef)
	}

	
	l <- list(param = param, mconc = mconc, dag = dag)
	return(l)
}

#' Sample a coefficient matrix of iid coefficients.
#'
#' @rdname iid
#'
#' @param p Positive integer, Number of dimension.
#' @param rmcoef function, the random number generator for the non-zero entries.
#' @param edges Matrix containing the edges of the Gaussian Bayesian network (from, to)
#'
#' @details Sample a sparse coefficient matrix for a Gaussian Bayesian network.
#'
#' @return  mcoef a lower triangular matrix with ones on the diagonal, and entry
#' (i, j) = 0 if (j, i) is not a row in `edges`.
#'
.rcoef_iid <- function(p = 100,
                       rmcoef = rnorm, 
					   edges)
{

	mcoef <- Matrix::sparseMatrix(i = edges[, 2], j = edges[, 1], x =
						  rmcoef(nrow(edges)), dims = c(p, p), triangular = TRUE)

	diag(mcoef) <- 1
	return(mcoef)
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
						 edges) 
{	

	mcoef <- sparseMatrix(i = edges[, 2], j = edges[, 1], dims = c(p, p),
						  triangular = TRUE, x = rep(1, nrow(edges)))
	theta <- drop0(tril(Matrix(pi/2, nrow = p, ncol = p), k = -1))

	diag(mcoef) <- 1

	for (icol in 1:(p - 1)) {
		for (irow in (icol + 1):p) {
			if (mcoef[irow, icol] != 0) {
				theta[irow, icol] <- .rsin(n = 1, k = p - icol, method = method)
				mcoef[irow, icol] <- cos(theta[irow, icol])
			} 
		}
		if (icol >= 2) {
			for (j in 1:(icol - 1)) {
				mcoef[icol:p, icol] <- mcoef[icol:p, icol] * sin(theta[icol:p, j])
			}
		}	
	}
	mcoef[p, p] <- prod(sin(theta[p, 1:(p - 1)]))
	l <- list(mcoef = mcoef, mangl = theta)
	return(l)
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


#'
#'@export
directUnifSampling <- function(dag, h=100, eps = 0.001){
  isCh<- is_chordal(dag,fillin=T)
  if (!isCh$chordal){
    warning("Can't sample uniformly for non chordal graph")
    dag <- add_edges(dag,edges = isCh$fillin)
  }
  edges <- as_edgelist(dag)
  v <- V(dag)
  p <- length(v)
  k <- rep(1,p) + degree(dag, mode = "in")
  U <- .rcoef_iid(p = p, edges = edges) 
  U <- t(U)
  U <- t(apply(U, MARGIN = 1, FUN = function(u) return(u/sqrt(sum(u^2))) ))
  J <- .ljac(diag(U), k)
  print(J)
  for (i in (1:h)){
    E <- t(.rcoef_iid(p = p, edges = edges, rmcoef = function(n){
      return(rnorm(n, mean = 0, sd = eps))
    }))
    diag(E)<-rnorm(p,mean=0,sd=eps)
    Up <- U + E
    Up <- t(apply(Up, MARGIN = 1, FUN = function(u) return(u/sqrt(sum(u^2))) ))
    Jp <- .ljac(diag(Up), k)
    print(Jp-J)
    if (log(runif(1)) <= (Jp - J) ){
      U <- Up
      J <- Jp
      print("step")
    }
  }
  return(U)
}


#' without checking
.ljac <- function(u, k){
  return(sum(k*log(abs(u))))
}



#' sampling on sphere proportionally to a power of the first coordinate
#' 
#' @param N sample size
#' @param n dimension
#' @param k exponent of the density
#' @param h heating phase size
#' @param s step
#' 
#'  Metropolis-Hasting algorithm to sample in the n dimensional semi-sphere (x_1>0)  
#' @export
.sphereSample <- function(N = 1, n, k = 1, h = 100, s = 1, eps = 0.01,returnAll=FALSE){
  Tot <- h + (N-1)*s + 1 #total number of iteration of MH
  Sample <- matrix(nrow = Tot, ncol = n) #obj initialization 
  Sample[1,] <- rnorm(n = n, mean = 0, sd = 1) #first point
  Sample[1,1]<-abs(Sample[1,1]) #absolute value first component (has to be positive)
  Sample[1,] <- Sample[1,]/sqrt(sum(Sample[1, ] ^ 2)) #normalization 
 for (i in 2:Tot){
    prop <- Sample[i-1,] + rnorm(n = n, mean = 0, sd = eps) # perturbate previous sample
    prop <- prop / sqrt(sum(prop ^ 2)) #normalize propozed
    #if (runif(1) <= sign(prop[1])*abs(prop[1])^k / Sample[i - 1 , 1]^k){
    if ((prop[1]>0) && (log(runif(1)) <= k*log((prop[1])) - k*log(Sample[i - 1 , 1]))){
      Sample[i,] <- prop
    }else{
      Sample[i,] <- Sample[i - 1 , ]
    }
  }
  if (!returnAll){
    Sample <- Sample[seq(from = h+1, to = Tot, by = s),]
  }
  
  return(Sample)
}


#'
#' @export
.sampleU <- function(N = 1, dag, h = 100, s = 1, eps = 0.1){
  #isCh<- is_chordal(dag,fillin=T)
  #if (!isCh$chordal){
    #warning("Can't sample uniformly for non chordal graph")
    #return(NULL)
  #}
  p <- length(V(dag)) 
  edges <- as_edgelist(dag)
  U <- array(dim=c(p,p,N),data = 0)
  
  u <-  t(as.matrix(Matrix::sparseMatrix(i = edges[, 2], j = edges[, 1], dims = c(p, p),
                             triangular = TRUE, x = rep(1, nrow(edges)))))
  diag(u)<-1
  U[,,1] <- u
  U[p,p,1:N] <- 1
  ch <- degree(dag,mode = "out")
  pa <- degree(dag, mode="in")
  for (i in 1:(p-1)){
    su <- .sphereSample(N = N,n = ch[i]+1 ,k = pa[i]+1,h = h, s = s,eps=eps)
    U[i,U[i,,1]>0,1:N] <- t(su)
  }
  return(U)
}


#' rUnifChordalConc
#' 
#' @param  N 
#' @param dag directed chordal acyclic graph, use the igraph package graph class or alternatively 
#'            a matrix that will be interpreted as an adjacency matrix. 
#' @param return.minvector logical, if TRUE the minimimal vector representation is returned (useful to plot in the elliptope)
#' @param ... additional parameters
#' 
#' @export
 rUnifChordalConc <- function(N=1, dag=NULL, return.minvector = FALSE, ... ){
   if (is.matrix(dag)){
     dag <- igraph::graph_from_adjacency_matrix(adjmatrix = dag, mode = directed)
   }
   isCh<- is_chordal(dag,fillin=T)
   if (!isCh$chordal){
     warning("Can't sample uniformly for non chordal graph")
     #return(NULL)
   }

   sU <- .sampleU(N = N, dag = dag, ... )
   vsC <- apply(sU,MARGIN = 3, function(U) return(U%*%t(U)) )
   sC <- array(data=vsC, dim=dim(sU))
   if (return.minvector){
     mv <- apply(sC, MARGIN = 3, function(m){
       return(m[upper.tri(m)])
     })
     return(t(mv))
   }else{
     return(sC)     
   }

 }


