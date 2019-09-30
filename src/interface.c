#include "interface.h"
#include "port.h"
#include "error.h"

SEXP C_port (SEXP R_madj, SEXP R_Q) {
	int *dims = NULL;
	unsigned int i = 0, N = 0, m_dim = 0, p = 0;
	SEXP R_mort = NULL, R_res = NULL, R_dims = NULL;
	double *mort, *madj, *Q, *res;
	gmat_errno_t gmat_errno;

	/* Initialize problem dimensions */
	/* Arrays always have dim attribute so no null check needed */
	R_dims = getAttrib(R_Q, R_DimSymbol);
	if (length(R_dims) != 3) {
		error("Expected 3 dimensions for Q array.\n");
	}
	dims = INTEGER(R_dims);
	if (dims[0] != dims[1]) {
		error("Q factors must be square.\n");
	}
	p = dims[0]; N = dims[2]; m_dim = p * p;

	/* Initialize R objects */
	R_mort = PROTECT(alloc3DArray(REALSXP, p, p, N));
	R_res = PROTECT(alloc3DArray(REALSXP, p, p, N));
	mort = REAL(R_mort);
	res = REAL(R_res);

	/* Get argument pointers */
	madj = REAL(R_madj); Q = REAL(R_Q);

	for (i = 0; i < N; i++) {
		gmat_errno = gram_schmidt_sel(mort + m_dim*i, madj, Q + m_dim*i, p);
		if (gmat_errno != GMAT_OK) {
			UNPROTECT(2);
			error("%s.\n", gmat_strerror(gmat_errno));
		}
		gmat_errno = crossproduct(res + m_dim*i, mort + m_dim*i, madj, p);
		if (gmat_errno != GMAT_OK) {
			UNPROTECT(2);
			error("%s\n", gmat_strerror(gmat_errno));
		}
	}

	UNPROTECT(2);
	return R_res;
}
