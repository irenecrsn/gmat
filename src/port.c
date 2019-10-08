#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "port.h"
#include "error.h"


static void port(double *mort, double *madj, double *Q, unsigned int p);

/*
 * Performs the crossproduct of two matrices in column major format.
 *
 * @param res Crossproduct matrix result (t(Q) * Q)
 * @param Q Matrix which will be multiplied
 * @param madj For hard-coding the zeros in the result (otherwise, extremely
 * small values are returned in such entries instead of zeros, because of the
 * poor numeric properties of Gram Schmidt orthogonalization).
 * @param p Matrix dimension (p x p)
 */
static void crossproduct (double *res, double *Q,
		double *madj, unsigned int p);

/*
 * Modified Gram Schmidt algorithm
 *
 * @param res Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param p Dimension of vectors (columns of span)
 */
static void gram_schmidt (double *res, double **span,
		unsigned int nvec, unsigned int p);

/*
 * Selective Gram Schmidt algorithm
 */
int port_sample(double *res, double *madj,	double *Q, unsigned int p,
	unsigned int N) {

	double mort[p * p];

	if (res == NULL || madj == NULL || Q == NULL) {
		return PORT_ENULL;
	}

	for (int n = 0; n < N; n++) {
		/* Partial orthogonalization of the initial Q factors */
		port(mort, madj, Q + n*p*p, p);
		/* Crossproduct t(Q) * Q of the resulting factors */
		crossproduct(res + n*p*p, mort, madj, p);
	}

	return PORT_OK;
}

static void port(double *res, double *madj, double *Q, unsigned int p) {

	double *span_sel[p], ort_base[p * p];
	unsigned int n_span = 0, jp = 0;

	memcpy(res, Q, sizeof(double) * p * p);

	for (int j = 0; j < p; j++) {

		jp = j * p;
		n_span = 0;

		for (int i = 0; i < j; i++) {
			if (madj[jp + i] == 0) {
				span_sel[n_span] = res + i * p;
				n_span++;
			}
		}
		span_sel[n_span] = res + jp;
		n_span++;

		/* we orthonormalize the span obtained for the current row */
		gram_schmidt(ort_base, span_sel, n_span, p);
		/* the last vector of the orthonormalized span is the new j-th column */
		memcpy(res + jp, ort_base + (n_span - 1)*p, sizeof(double) * p);
	}
}

/*
 * Modified Gram Schmidt algorithm
 */
static void gram_schmidt (double *res, double **span,
		unsigned int nvec, unsigned int p) {

	double norm = 0, dot_prod = 0;
	unsigned int ip = 0, jp = 0;

	for (int i = 0; i < nvec; i++) {
		ip = i * p;
		memcpy(res + ip, span[i], sizeof(double) * p);
		for (int j = 0; j < i; j++) {
			jp = j * p;
			/* Orthogonal projection */
			dot_prod = 0;
			for (int k = 0; k < p; k++) {
				dot_prod += (res[ip + k] * res[jp + k]);
			}
			for (int k = 0; k < p; k++) {
				res[ip + k] -= (dot_prod * res[jp + k]);
			}
		}
		/* Normalization */
		norm = 0;
		for (int k = 0; k < p; k++) {
			norm += (res[ip + k] * res[ip + k]);
		}
		norm = 1 / sqrt(norm);
		for (int k = 0; k < p; k++) {
			res[ip + k] *= norm;
		}
	}
}

/*
 * Performs the crossproduct of two matrices in column major format.
 */
static void crossproduct (double *res, double *mort, double *madj,
							 unsigned int p) {
	unsigned int jp = 0, ip = 0;
	double sum = 0;

	/* Upper triangle first */
	for (int j = 0; j < p; j++) {
		jp = j * p;
		for (int i = 0; i < j; i++) {
			if (madj[jp + i] == 0) {
				res[jp + i] = 0; /* Hard-code 0s for missing edges */
			} else { /* Crossproduct */
				sum = 0;
				ip = i * p;
				for (int k = 0; k < p; k++) {
					sum += (mort[jp + k]*mort[k + ip]);
				}
				res[jp + i] = sum;
			}
		}
	}

	/* Lower triangle == transpose of upper triangle (crossproduct) */
	for (int j = 0; j < p; j++) {
		res[j * p + j] = 1; /* Diagonal = 1 (normalized vectors) */
		for (int i = j + 1; i < p; i++) {
			res[j * p + i] = res[i * p + j];
		}
	}
}
