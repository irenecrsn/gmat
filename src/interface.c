#include "interface.h"
#include "gmn_sampling.h"

SEXP C_port (SEXP madj, SEXP Q) {
	int *dims;
	unsigned int i = 0, N, matrix_dim;;
	SEXP mort = NULL, R_dims = NULL;
	double *p_mort, *p_madj, *p_Q;;

	R_dims = getAttrib(Q, R_DimSymbol);
	dims = INTEGER(R_dims);

	N = dims[2];
	mort = PROTECT(alloc3DArray(REALSXP, dims[0], dims[1], N));
	p_mort = REAL(mort);
	p_madj = REAL(madj);
	p_Q = REAL(Q);
	matrix_dim = dims[0] * dims[1];

	for (i = 0; i < N; i++) {
		gram_schmidt_sel(p_mort + matrix_dim*i, p_madj, p_Q + matrix_dim*i, dims[0]);
		crossprod(p_mort + matrix_dim*i,
	}

	UNPROTECT(1);
	return mort;
}
