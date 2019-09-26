#include "interface.h"
#include "gmn_sampling.h"
#include "utils.h"

SEXP C_port (SEXP madj, SEXP Q) {
	int *dims;
	unsigned int i = 0, N, matrix_dim;;
	SEXP mort = NULL, R_dims = NULL, res = NULL;
	double *p_mort, *p_madj, *p_Q, *p_res;

	R_dims = getAttrib(Q, R_DimSymbol);
	dims = INTEGER(R_dims);

	N = dims[2];
	mort = PROTECT(alloc3DArray(REALSXP, dims[0], dims[1], N));
	res = PROTECT(alloc3DArray(REALSXP, dims[0], dims[1], N));
	p_mort = REAL(mort);
	p_res = REAL(res);
	p_madj = REAL(madj);
	p_Q = REAL(Q);
	matrix_dim = dims[0] * dims[1];

	for (i = 0; i < N; i++) {
		gram_schmidt_sel(p_mort + matrix_dim*i, p_madj, p_Q + matrix_dim*i, dims[0]);
		crossprod(p_mort + matrix_dim*i, p_res + matrix_dim*i, dims[0],
		dims[1]);
	}

	UNPROTECT(2);
	return res;
}
