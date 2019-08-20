
### implemented with QR decomposition 
.rgmn_sqrt_qr <- function(rmsqrt, ug) {

        p <- length(V(ug))

	msqrt <- matrix(nrow = p, ncol = p, data = rmsqrt(p^2))

	msqrt <- .ort_selective_qr(ug, msqrt)
	
	return(msqrt %*% t(msqrt))
}

.ort_selective_qr <- function(ug, mort) {

       madj <- igraph::as_adjacency_matrix(ug, type = "both", sparse = FALSE)
       p <- nrow(madj)
       degrees <- degree(ug)
       order <- order(degrees, decreasing = FALSE)
      # order <- 1:p
        n_zeros <- length(degrees[degrees == 0])  
       
       ### first we orthogonalize all the rows of the disconnected nodes (+1)
        if (n_zeros > 0 ){
          qrdec <- qr(t(mort[order[1:(n_zeros+1)],]))
          mort[order[1:(n_zeros+1)],] <- t(qr.Q(qrdec)) 
        }
       
        ## now we orthogonalize the other rows following the order
	for (i in (n_zeros+2):p) {
	
                row_iort <- which(madj[order[i], ] == FALSE)

		row_iort <- row_iort[ row_iort %in% order[1:(i-1)]]
		row_nort <- length(row_iort)
		if (row_nort > 1) {
		    # mort[order[i], ] <- mort[order[i], ] -
                    # t( mort[row_iort,] )  %*% (solve( mort[row_iort,]%*% t(
                    # mort[row_iort,]) ) %*%  mort[row_iort,]  %*% mort[order[i],])

                    qrdec <- qr(t(mort[c(row_iort,order[i]),]))
                    mort[order[i],] <- t(qr.Q(qrdec))[row_nort + 1,]

                   # dec <- qr(t(mort[row_iort,]))
                   # Rr <- qr.R(dec)
                   #  mort[order[i], ] <- mort[order[i], ] -
                   #  t( mort[row_iort,] )  %*% (solve( t(Rr)%*% (Rr )) %*%  mort[row_iort,]  %*% mort[order[i],])
                  
                    
		} else if (row_nort == 1) {
                    mort[order[i],] <- mort[order[i],] - mort[row_iort,]* sum(mort[order[i],]* mort[row_iort,])/ sum(mort[row_iort,]^2)
		}
		}

	return (mort)
}

