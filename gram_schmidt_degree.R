

require(Matrix)
require(igraph)

.rgmn_sqrt_new <- function(rmsqrt, ug) {

        p <- length(V(ug))

	msqrt <- matrix(nrow = p, ncol = p, data = rmsqrt(p^2))

	msqrt <- .ort_selective(ug, msqrt)
	
	return(msqrt %*% t(msqrt))
}



.ort_selective <- function(ug, mort) {
       
       madj <- igraph::as_adjacency_matrix(ug, type = "both", sparse = FALSE)
       p <- nrow(madj)
       degrees <- degree(ug)
       order <- order(degrees, decreasing = FALSE)
       n_zeros <- length(degrees[degrees == 0]) + 1 
      
        if (n_zeros > 1){
          mort[order[1:n_zeros],] <- .gram_schmidt(mort[order[1:n_zeros],])
        # qrdec <- qr(t(mort[order[1:(n_zeros+1)],]))
        # mort[order[1:(n_zeros+1)],] <- t(qr.Q(qrdec))
        }

	for (i in (n_zeros+1):p) {

		row_iort <- which(madj[order[i], ] == FALSE)
		row_iort <- row_iort[ row_iort %in% order[1:(i-1)]]
		row_nort <- length(row_iort)

		if (row_nort > 0) {
			row_ort_base <- .gram_schmidt(mort[row_iort, ],k=
                        max(n_zeros-1,1))
			if (row_nort == 1) {
				mort[i, ] <- mort[i, ] - .proj(mort[i, ], row_ort_base)
			} else {
				for (j in 1:row_nort) {
					mort[i, ] <- mort[i, ] - .proj(mort[i, ], row_ort_base[j, ])
				}
			}
		}
	}

	return (mort)
}

.gram_schmidt <- function(span,k=1) {

	if (class(span) == "numeric") {
		return (span)
	}

	for (i in (k+1):nrow(span)) {
		for (j in 1:(i - 1)) {
			span[i, ] <- span[i, ] - .proj(span[i, ], span[j, ])
		}
	}

	return (span)
}

# orthogonal projection of v on u
.proj <- function(v, u) {
	dot_uv <- sum(u * v)
	dot_uu <- sum(u * u)
	
	return ((dot_uv/dot_uu)*u)
}
