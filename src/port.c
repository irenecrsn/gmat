#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef NDEBUG
#define NDEBUG
#endif
#include <assert.h>

#include "port.h"
#include "error.h"

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
 * Gram Schmidt algorithm
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param p Dimension of vectors (columns of span)
 */
static void gram_schmidt (double *span_ort, double **span,
		unsigned int nvec, unsigned int p);

static void proj_ort (double *v_proj_u, double *v,
		double *u, unsigned int p);

/*
 * Selective Gram Schmidt algorithm
 */
int port(double *res, double *madj,	double *Q, unsigned int p) {

	double *span_sel[p], *ort_base = NULL, *mort = NULL;
	unsigned int i = 0, j = 0;
	unsigned int n_span = 0, jp = 0, m_dim = 0;

	if (res == NULL || madj == NULL || Q == NULL) {
		return GMAT_ENULL;
	}

	m_dim = p * p;
	ort_base = calloc(m_dim, sizeof(double));
	mort = calloc(m_dim, sizeof(double));

	if (ort_base == NULL || mort == NULL) {
		free(ort_base); ort_base = NULL;
		free(mort); mort = NULL;
		return GMAT_ENOMEM;
	}
	memcpy(mort, Q, sizeof(double) * m_dim);

	/* Partial orthogonalization of the initial Q factors */
	for (j = 0; j < p; j++) {

		jp = j * p;
		n_span = 0;

		for (i = 0; i < j; i++) {
			if (madj[jp + i] == 0) {
				span_sel[n_span] = mort + i * p;
				n_span++;
			}
		}
		span_sel[n_span] = mort + jp;
		n_span++;

		/* we orthonormalize the span obtained for the current row */
		gram_schmidt(ort_base, span_sel, n_span, p);
		/* the last vector of the orthonormalized span is the new j-th column */
		memcpy(mort + jp, ort_base + (n_span - 1)*p, sizeof(double) * p);
	}
	/* Crossproduct t(Q) * Q of the resulting factors */
	crossproduct(res, mort, madj, p);

	free(ort_base); ort_base = NULL;
	free(mort); mort = NULL;
	return GMAT_OK;
}


/*
 * Performs the crossproduct of two matrices in column major format.
 */
static void crossproduct (double *res, double *mort, double *madj,
							 unsigned int p) {
	unsigned int i = 0, j = 0, k = 0, jp = 0, ip = 0;
	double sum = 0;

	assert(mort != NULL); assert(madj != NULL);

	/* Upper triangle first */
	for (j = 0; j < p; j++) {
		jp = j * p;
		for (i = 0; i < j; i++) {
			if (madj[jp + i] == 0) {
				res[jp + i] = 0; /* Hard-code 0s for missing edges */
			} else { /* Crossproduct */
				sum = 0;
				ip = i * p;
				for (k = 0; k < p; k++) {
					sum += (mort[jp + k]*mort[k + ip]);
				}
				res[jp + i] = sum;
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
}

/*
 * Gram Schmidt algorithm
 */
static void gram_schmidt (double *span_ort, double **span,
		unsigned int nvec, unsigned int p)
{
	double v_proj[p];
	unsigned int i = 0, j = 0, k = 0, ip = 0;
	double norm = 0;

	assert(span_ort != NULL); assert(span != NULL); assert(nvec != NULL);

	for (i = 0; i < nvec; i++) {
		memcpy(span_ort + i*p, span[i], sizeof(double) * p);
		ip = i * p;
		for (j = 0; j < i; j++) {
			proj_ort(v_proj, span_ort + ip, span_ort + j*p, p);
			for (k = 0; k < p; k++) {
				span_ort[ip + k] -= v_proj[k];
			}
		}
		/* we normalize the resulting vector */
		norm = 0;
		for (k = 0; k < p; k++) {
			norm += span_ort[ip + k] * span_ort[ip + k];
		}
		norm = 1 / sqrt(norm);
		for (k = 0; k < p; k++) {
			span_ort[ip + k] *= norm;
		}
	}
}

/*
 * Orthogonal projection of v onto direction u.
 * Vector u is assumed to already be normalized.
 */
static void proj_ort (double *v_proj_u, double *v, double *u, unsigned int p)
{
	unsigned int i = 0;
	double dot_uv = 0;

	assert(u != NULL); assert(v != NULL); assert(v_proj_u != NULL);

	for (i = 0; i < p; i++) {
		dot_uv += (u[i] * v[i]);
	}

	for (i = 0; i < p; i++) {
		v_proj_u[i] = dot_uv * u[i];
	}
}
