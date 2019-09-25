#include "interface.h"
#include "gmn_sampling.h"

SEXP C_gram_schmidt_sel (SEXP madj, SEXP mcov) {
	int n_col = ncols(madj);
	SEXP mort;
	double *p_mort;

	mort = PROTECT(allocMatrix(REALSXP, n_col, n_col));
	p_mort = REAL(mort);
    
	gram_schmidt_sel(p_mort, REAL(madj), REAL(mcov), n_col);
    
	UNPROTECT(1);
	return mort;
}
