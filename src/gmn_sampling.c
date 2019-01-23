#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gmn_sampling.h"


int compare ( const void *pa, const void *pb )
{
  const int *a = pa;
  const int *b = pb;
  if(a[0] < b[0]){
    return -1;
  }else if (a[0] > b[0]){
    return 1;
  }else{
    if(a[1] < b[1]){
      return -1;
    }else if (a[1] > b[1]){
      return 1;
    }else{
      if(a[2] < b[2]){
        return -1;
      }else if (a[2] > b[2]){
        return 1;
      }else{
        if(a[3] < b[3]){
          return -1;
        }else if (a[3] > b[3]){
          return 1;
        }else{
          return 0; 
        }
      }
    } 
  }
}

/*
 * Selective Gram Schmidt algorithm 
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param dim Dimension of vectors (columns of span)
 */
int gram_schmidt_sel (double *mort, int *madj, double *mcov, 
                      unsigned int *dim) {
  double **span_sel = NULL, **ort_base = NULL;
  double *v_proj = NULL;
  double nn = 0;
  int allright = 0;
  int maps[dim[0]][4], ix[dim[0] + 1], cc[dim[0]];
  unsigned int i = 0, j = 0, k = 0, skip = 0;
  unsigned int n_span = 0, i_current = 0, nzeros = 0;
  ix[0] = -1;
  
  if (mort == NULL || madj == NULL || mcov == NULL || dim == NULL) {
    return -1;
  }
  
  
  if ((v_proj = calloc(dim[0], sizeof(double))) == NULL) {
    return -1;
  }
  
  if ((span_sel = calloc(dim[0], sizeof(double *))) == NULL) {
    free(v_proj); v_proj = NULL;
    return -1;
  }
  
  if ((ort_base = calloc(dim[0], sizeof(double *))) == NULL) {
    free(span_sel); span_sel = NULL;
    free(v_proj); v_proj = NULL;
    return -1;
  }
  
  for (i = 0; i < dim[0]; i++) {
    ort_base[i] = NULL;
  }
  
  for (i = 0; i < dim[0]; i++) {
    if ((ort_base[i] = calloc(dim[0], sizeof(double))) == NULL) {
      for (j = 0; j < i; j++) {
        free(ort_base[j]); ort_base[j] = NULL;
      }
      free(v_proj); v_proj = NULL;
      free(span_sel); span_sel = NULL;
      free(ort_base); ort_base = NULL;
      return -1;
    }
  }
  
  for (i = 0; i < dim[0]; i++){
    maps[i][1] = -1;
  }
  
  k = 0;
  /* compute degrees, connected components and size of cc and store the index*/
  for (i = 0; i < dim[0]; i++) {
    if (maps[i][1] < 0){ /*new connected components*/
      k++; /* increment the connected components*/
      maps[i][1] = k; 
    } 
    cc[maps[i][1]]++;
    maps[i][2] = 0;
    maps[i][3] = i; /* store the index*/
  for (j = 0; j < dim[0]; j++) {
    if ( (madj[i * dim[0] + j] > 0) && (j != i) ){
      maps[i][2]++;  /* increase count of degree*/             
      maps[j][1] = maps[i][1];  /*propagate the index of the conn comp */
    }
  }
  }
  /*copy size of con comps*/
  for (i = 0; i < dim[0]; i++){ 
    maps[i][0] = cc[maps[i][1]];
  }


  /* sort the maps */
  qsort(maps, dim[0], sizeof(int) * 4, compare);
  
  nzeros = 0;
  while (maps[nzeros][2] == 0){
    nzeros++;
  }
  
  /* we ortogonalize the disconnected ones */
  n_span = 0;
  for (j = 0; j < nzeros; j++) {
    span_sel[n_span] = mcov + maps[j][3] * dim[0];
    n_span++;
  }
  gram_schmidt(ort_base, span_sel, &n_span, dim, 0);
  /* and we copy them in the result */
  for (j = 0; j < nzeros; j++) {
    for (k = 0; k < dim[0]; k++) {
      mort[maps[j][3] * dim[0] + k] = ort_base[j][k];
    }
  }
  
  /* from now on in the first nzeros components of ort_base there
   * are orthogonal vectors for the disconnected nodes
   */
  
  /* now the remaining */
  for (i = nzeros; i < dim[0]; i++) {
    i_current = maps[i][3] * dim[0];
    memcpy(mort + i_current, mcov + i_current, sizeof(double) * dim[0]);
    n_span = nzeros;
    skip = nzeros;
    allright = 1;
    for (j = nzeros; j < i; j++) {
      if (madj[i_current + maps[j][3]] == 0) {
        if ( (ix[n_span - nzeros] ==  maps[j][3]) && (allright == 1) ){
          skip++;
        }else{
          ix[n_span - nzeros] = maps[j][3];
          allright = 0;
        }
        span_sel[n_span] = mort + maps[j][3] * dim[0];
        n_span++;
      }
    }
    ix[n_span - nzeros] = maps[i][3];
    ix[n_span - nzeros + 1] = -1;
    span_sel[n_span] = mort + i_current;
    n_span++;
    gram_schmidt(ort_base, span_sel, &n_span, dim, skip);
    for (k = 0; k < dim[0]; k++) {
      mort[i_current + k] = ort_base[n_span - 1][k];
    }
  }
  
  free(v_proj); v_proj = NULL;
  free(span_sel); span_sel = NULL;
  for (i = 0; i < dim[0]; i++) {
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
 * @param skip skip the first skip vector of span (already orthogonal)
 */
int gram_schmidt (double **span_ort, double **span, 
                  unsigned int *nvec, unsigned int *dim, unsigned int skip) 
{
  double *v_proj = NULL;
  unsigned int i = 0, j = 0, k = 0;
  double nn = 0;
  
  if (span_ort == NULL || span == NULL || nvec == NULL || dim == NULL) {
    return -1;
  }
  
  for (i = skip; i < nvec[0]; i++) {
    memcpy(span_ort[i], span[i], sizeof(double) * dim[0]);
  }
  
  if ((v_proj = calloc(dim[0], sizeof(double))) == NULL) {
    return -1;
  }
  
  
  for (i = skip; i < nvec[0]; i++) {
    for (j = 0; j < i; j++) {
      proj_ort(v_proj, span_ort[i], span_ort[j], dim);
      for (k = 0; k < dim[0]; k++) {
        span_ort[i][k] -= v_proj[k];
      }
    }
    nn = 0;
    for (k = 0; k < dim[0]; k++) {
      nn += span_ort[i][k] * span_ort[i][k] ;
    }
    nn = 1 / sqrt(nn);
    for (k = 0; k < dim[0]; k++) {
      span_ort[i][k] = span_ort[i][k] * nn ;
    }
  }
  
  free(v_proj); v_proj = NULL;
  
  return 0;
}



/* 
 * Orthogonal projection of v onto u, the vector u is assumed to be normalized
 */
int proj_ort (double *v_proj_u, double *v, double *u, unsigned int *dim)
{
  unsigned int i = 0;
  double dot_uv = 0, dot_uu = 0;
  
  if (v_proj_u == NULL || v == NULL || u == NULL || dim == NULL) {
    return -1;
  }
  
  for (i = 0; i < dim[0]; i++) {
    dot_uv += (u[i] * v[i]);
  }
  
  for (i = 0; i < dim[0]; i++) {
    v_proj_u[i] = dot_uv * u[i];
  }
  
  return 0;
}

