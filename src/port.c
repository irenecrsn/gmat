#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NDEBUG
#include <assert.h>

#include "port.h"
#include "error.h"

static int gram_schmidt (double *span_ort, double **span,
		const unsigned int nvec, const unsigned int p);

static void proj_ort (double *v_proj_u, const double *v,
		const double *u, const unsigned int p);

/*
 * Selective Gram Schmidt algorithm
 */
int gram_schmidt_sel (double *mort, const double *madj,
		const double *mcov, const unsigned int p) {

	double *v_proj = NULL, *ort_base = NULL, **span_sel = NULL;
	unsigned int i = 0, j = 0;
	unsigned int n_span = 0, j_current = 0;

	if (mort == NULL || madj == NULL || mcov == NULL) {
		return GMAT_ENULL;
	}

	v_proj = calloc(p, sizeof(double));
	span_sel = calloc(p, sizeof(double *));
	ort_base = calloc(p * p, sizeof(double));

	if (v_proj == NULL || span_sel == NULL || ort_base == NULL) {
		free(v_proj); v_proj = NULL;
		free(span_sel); span_sel = NULL;
		free(ort_base); ort_base = NULL;
		return GMAT_ENOMEM;
	}

	for (j = 0; j < p; j++) {

		j_current = j * p;
		memcpy(mort + j_current, mcov + j_current, sizeof(double) * p);
		n_span = 0;

		for (i = 0; i < j; i++) {
			if (madj[j_current + i] == 0) {
				span_sel[n_span] = mort + i * p;
				n_span++;
			}
		}
		span_sel[n_span] = mort + j_current;
		n_span++;

		/* we orthonormalize the span obtained for the current row */
		gram_schmidt(ort_base, span_sel, n_span, p);
		for (i = 0; i < p; i++) {
			mort[j_current + i] = ort_base[(n_span - 1)*p + i];
		}
	}

	free(v_proj); v_proj = NULL;
	free(span_sel); span_sel = NULL;
	free(ort_base); ort_base = NULL;
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

/*
 * Gram Schmidt algorithm
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param p Dimension of vectors (columns of span)
 */
static int gram_schmidt (double *span_ort, double **span,
		const unsigned int nvec, const unsigned int p)
{
	double *v_proj = NULL;
	unsigned int i = 0, j = 0, k = 0, i_current = 0;
	double norm = 0;

	assert(span_ort != NULL); assert(span != NULL); assert(nvec != NULL);

	for (i = 0; i < nvec; i++) {
		memcpy(span_ort + i*p, span[i], sizeof(double) * p);
	}

	if ((v_proj = calloc(p, sizeof(double))) == NULL) {
		return GMAT_ENOMEM;
	}

	for (i = 0; i < nvec; i++) {
		i_current = i * p;
		for (j = 0; j < i; j++) {
			proj_ort(v_proj, span_ort + i_current, span_ort + j*p, p);
			for (k = 0; k < p; k++) {
				span_ort[i_current + k] -= v_proj[k];
			}
		}
		/* we normalize the resulting vector */
		norm = 0;
		for (k = 0; k < p; k++) {
			norm += span_ort[i_current + k] * span_ort[i_current + k];
		}
		norm = 1 / sqrt(norm);
		for (k = 0; k < p; k++) {
			span_ort[i_current + k] = span_ort[i_current + k] * norm;
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
		const unsigned int p)
{
	unsigned int i = 0;
	double dot_uv = 0;

	assert(u != NULL); assert(v != NULL);
	assert(v_proj_u != NULL);

	for (i = 0; i < p; i++) {
		dot_uv += (u[i] * v[i]);
	}

	for (i = 0; i < p; i++) {
		v_proj_u[i] = dot_uv * u[i];
	}
}
