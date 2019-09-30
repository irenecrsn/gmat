#ifndef PORT_H
#define PORT_H

/*
 * Partial orthogonalization sampling method
 *
 * @param mort Result matrix (in column major)
 * @param madj Adjacency matrix of the graph (gives the orthogonal pairs of
 * vectors)
 * @param mcov Matrix containing as columns the initial matrix sample
 * @param p Dimension of the matrices (p x p)
 */
int port(double *res, const double *madj,
		const double *mcov, const unsigned int p);

#endif
