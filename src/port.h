#ifndef PORT_H
#define PORT_H

/*
 * Selective Gram Schmidt algorithm
 *
 * @param mort Result matrix (in column major)
 * @param madj Adjacency matrix of the graph (gives the orthogonal pairs of
 * vectors)
 * @param mcov Matrix containing as columns the vectors to orthogonalize
 * @param p Dimension of the matrices (p x p) 
 */
int gram_schmidt_sel (double *mort, const double *madj,
		const double *mcov, const unsigned int p);

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
int crossproduct (double *res, const double *Q,
		const double *madj, const unsigned int p);

#endif
