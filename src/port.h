#ifndef PORT_H
#define PORT_H

/*
 * Selective Gram Schmidt algorithm
 *
 * @param span_ort Matrix with rows containing the orthogonal vectors
 * @param span Matrix with rows containing the vectors to orthogonalize
 * @param nvec Number of vectors to orthogonalize (rows of span)
 * @param dim Dimension of vectors (columns of span)
 */
int gram_schmidt_sel (double *mort, const double *madj,
		const double *mcov, const unsigned int dim);

int crossproduct (double *res, const double *mort,
		const double *madj, const unsigned int p);

#endif
