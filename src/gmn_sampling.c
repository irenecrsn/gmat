#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NDEBUG
#include <assert.h>

#include "gmn_sampling.h"

static int gram_schmidt (double **span_ort, double **span,
		const unsigned int nvec, const unsigned int dim);

static void proj_ort (double *v_proj_u, const double *v,
		const double *u, const unsigned int dim);

/*
 * Selective Gram Schmidt algorithm
 */
int gram_schmidt_sel (double *mort, const double *madj,
		const double *mcov, const unsigned int dim) {

	double **span_sel = NULL, **ort_base = NULL;
	double *v_proj = NULL;
	unsigned int i = 0, j = 0;
	unsigned int n_span = 0, j_current = 0;

	#define FREEALL \
	for (j = 0; j < i; j++) {\
		free(ort_base[j]); ort_base[j] = NULL;\
	}\
	free(v_proj); v_proj = NULL;\
	free(span_sel); span_sel = NULL;\
	free(ort_base); ort_base = NULL;

	if (mort == NULL || madj == NULL || mcov == NULL) {
		return GMAT_ENULL;
	}

	if ((v_proj = calloc(dim, sizeof(double))) == NULL) {
		return GMAT_ENOMEM;
	}

	if ((span_sel = calloc(dim, sizeof(double *))) == NULL) {
		FREEALL;
		return GMAT_ENOMEM;
	}

	if ((ort_base = calloc(dim, sizeof(double *))) == NULL) {
		FREEALL;
		return GMAT_ENOMEM;
	}

	for (i = 0; i < dim; i++) {
		ort_base[i] = NULL;
		if ((ort_base[i] = calloc(dim, sizeof(double))) == NULL) {
			FREEALL;
			return GMAT_ENOMEM;
		}
	}

	for (j = 0; j < dim; j++) {

		j_current = j * dim;
		memcpy(mort + j_current, mcov + j_current, sizeof(double) * dim);
		n_span = 0;

		for (i = 0; i < j; i++) {
			if (madj[j_current + i] == 0) {
				span_sel[n_span] = mort + i * dim;
				n_span++;
			}
		}
		span_sel[n_span] = mort + j_current;
		n_span++;

		/* we orthonormalize the span obtained for the current row */
		gram_schmidt(ort_base, span_sel, n_span, dim);
		for (i = 0; i < dim; i++) {
			mort[j_current + i] = ort_base[n_span - 1][i];
		}
	}

	FREEALL;
	return GMAT_OK;
}


int crossproduct (double * res, const double *mort, const double *madj,
							 const unsigned int p) {
	unsigned int i = 0, j = 0, k = 0, j_current = 0, i_current = 0;
	double sum = 0;

	if (mort == NULL || madj == NULL) {
		return GMAT_ENULL;
	}

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

	return GMAT_OK;
}

const char * gmat_strerror (const gmat_errno_t gmat_errno) {

	switch(gmat_errno) {
		case GMAT_ENULL:
			return "Unexpected NULL pointer";
		case GMAT_ENOMEM:
			return "Could not allocate more memory";
		default:
			return "Unknown error code";
	}
}

/*
 * Gram Schmidt algorithm
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param dim Dimension of vectors (columns of span)
 */
static int gram_schmidt (double **span_ort, double **span,
		const unsigned int nvec, const unsigned int dim)
{
	double *v_proj = NULL;
	unsigned int i = 0, j = 0, k = 0;
	double norm = 0;

	assert(span_ort != NULL); assert(span != NULL); assert(nvec != NULL);

	for (i = 0; i < nvec; i++) {
		memcpy(span_ort[i], span[i], sizeof(double) * dim);
	}

	if ((v_proj = calloc(dim, sizeof(double))) == NULL) {
		return GMAT_ENOMEM;
	}

	for (i = 0; i < nvec; i++) {
		for (j = 0; j < i; j++) {
			proj_ort(v_proj, span_ort[i], span_ort[j], dim);
			for (k = 0; k < dim; k++) {
				span_ort[i][k] -= v_proj[k];
			}
		}
		/* we normalize the resulting vector */
		norm = 0;
		for (k = 0; k < dim; k++) {
			norm += span_ort[i][k] * span_ort[i][k];
		}
		norm = 1 / sqrt(norm);
		for (k = 0; k < dim; k++) {
			span_ort[i][k] = span_ort[i][k] * norm;
		}
	}

	free(v_proj); v_proj = NULL;

	return GMAT_OK;
}

/*
 * Orthogonal projection of v onto direction u.
 * Vector u is assumed to already be normalized.
 */
static void proj_ort (double *v_proj_u, const double *v, const double *u,
		const unsigned int dim)
{
	unsigned int i = 0;
	double dot_uv = 0;

	assert(u != NULL); assert(v != NULL);
	assert(v_proj_u != NULL);

	for (i = 0; i < dim; i++) {
		dot_uv += (u[i] * v[i]);
	}

	for (i = 0; i < dim; i++) {
		v_proj_u[i] = dot_uv * u[i];
	}
}
