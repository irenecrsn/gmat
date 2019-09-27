#include "interface.h"
#include "gmn_sampling.h"

static int crossproduct (double *mort, double *madj, double *res,
							 unsigned int p);

SEXP C_port (SEXP R_madj, SEXP R_Q) {
	int *dims = NULL;
	unsigned int i = 0, N = 0, m_dim = 0, p = 0;
	SEXP R_mort = NULL, R_res = NULL;
	double *mort, *madj, *Q, *res;

	/* Initialize problem dimensions */
	/* Arrays always have dim attribute so no null check needed */
	dims = INTEGER(getAttrib(R_Q, R_DimSymbol));
	if (length(dims) != 3) {
		error("Expected 3 dimensions for Q array\n");
	}
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
		gram_schmidt_sel(mort + m_dim*i, madj, Q + m_dim*i, p);
		crossproduct(mort + m_dim*i, madj, res + m_dim*i, p);
	}

	UNPROTECT(2);
	return R_res;
}

static int crossproduct (double *mort, double *madj, double *res,
							 unsigned int p) {
	unsigned int i = 0, j = 0, k = 0, j_current = 0, i_current = 0;
	double sum = 0;

	/* Upper triangle first */
	for (j = 0; j < p; j++) {
		j_current = j * p;
		for (i = 0; i < j; i++) {
			if (madj[j_current + i] == 0) {
				res[j_current + i] = 0; /* Hard-code 0s for missing edges */
			} else { /* Crossproduct */
				sum = 0;
				i_current = i * p;
				for (k = 0; k < p; k++) {
					sum += (mort[j_current + k]*mort[k + i_current]);
				}
				res[j_current + i] = sum;
			}
		}
	}

	/* Lower triangle == transpose of upper triangle (crossproduct) */
	for (j = 0; j < p; j++) {
		res[j * p + j] = 1; /* Diagonal = 1 (normalized vectors) */
		for (i = j + 1; i < p; i++) {
			res[j * p + i] = res[i * p + j];
		}
	}

	return 0;
}
