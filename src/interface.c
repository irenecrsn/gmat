#include "interface.h"
#include "port.h"
#include "error.h"

SEXP C_port (SEXP R_madj, SEXP R_Q) {
	int *dims = NULL;
	unsigned int N = 0, p = 0;
	SEXP R_res = NULL, R_dims = NULL;
	double *madj, *Q, *res;
	port_errno_t port_errno;

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
	p = dims[0]; N = dims[2];

	/* Initialize R objects */
	R_res = PROTECT(alloc3DArray(REALSXP, p, p, N));
	res = REAL(R_res);

	/* Get argument pointers */
	madj = REAL(R_madj); Q = REAL(R_Q);

	port_errno = port_sample(res, madj, Q, p, N);
	if (port_errno != PORT_OK) {
		UNPROTECT(1);
		error("%s.\n", port_strerror(port_errno));
	}

	UNPROTECT(1);
	return R_res;
}
