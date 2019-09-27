#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NDEBUG
#include <assert.h>

#include "gmn_sampling.h"


/*
 * Selective Gram Schmidt algorithm 
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param dim Dimension of vectors (columns of span)
 */
int gram_schmidt_sel (double *mort, double *madj, double *mcov, 
		unsigned int dim) {
	double **span_sel = NULL, **ort_base = NULL;
	double *v_proj = NULL;
	unsigned int i = 0, j = 0;
	unsigned int n_span = 0, j_current = 0;
	
	if (mort == NULL || madj == NULL || mcov == NULL) {
		return -1;
	}

	if ((v_proj = calloc(dim, sizeof(double))) == NULL) {
		return -1;
	}

	if ((span_sel = calloc(dim, sizeof(double *))) == NULL) {
		free(v_proj); v_proj = NULL;
		return -1;
	}

	if ((ort_base = calloc(dim, sizeof(double *))) == NULL) {
		free(span_sel); span_sel = NULL;
		free(v_proj); v_proj = NULL;
		return -1;
	}

	for (i = 0; i < dim; i++) {
		ort_base[i] = NULL;
	}

	for (i = 0; i < dim; i++) {
		if ((ort_base[i] = calloc(dim, sizeof(double))) == NULL) {
			for (j = 0; j < i; j++) {
				free(ort_base[j]); ort_base[j] = NULL;
			}
			free(v_proj); v_proj = NULL;
			free(span_sel); span_sel = NULL;
			free(ort_base); ort_base = NULL;
			return -1;
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
		gram_schmidt(ort_base, span_sel, &n_span, dim);
		for (i = 0; i < dim; i++) {
			mort[j_current + i] = ort_base[n_span - 1][i];
		}
	}
	
	free(v_proj); v_proj = NULL;
	free(span_sel); span_sel = NULL;
	for (i = 0; i < dim; i++) {
		free(ort_base[i]); ort_base[i] = NULL;
	}
	free(ort_base); ort_base = NULL;

	return 0;
}


/*
 * Gram Schmidt algorithm 
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param dim Dimension of vectors (columns of span)
 */
int gram_schmidt (double **span_ort, double **span, 
		unsigned int *nvec, unsigned int dim) 
{
	double *v_proj = NULL;
	unsigned int i = 0, j = 0, k = 0;
	double norm = 0;

	assert(span_ort != NULL); assert(span != NULL); assert(nvec != NULL);

	for (i = 0; i < nvec[0]; i++) {
		memcpy(span_ort[i], span[i], sizeof(double) * dim);
	}

	if ((v_proj = calloc(dim, sizeof(double))) == NULL) {
		return -1;
	}

	for (i = 0; i < nvec[0]; i++) {
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

	return 0;
}

/* 
 * Orthogonal projection of v onto direction u.
 * Vector u is assumed to already be normalized.
 */
int proj_ort (double *v_proj_u, double *v, double *u, unsigned int dim)
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

	return 0;
}

