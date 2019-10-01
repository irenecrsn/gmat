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
 * @param N Sample size
 */
int port_sample(double *res, double *madj,	double *mcov, unsigned int p,
	unsigned int N);

#endif
