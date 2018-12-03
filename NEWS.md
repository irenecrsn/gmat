# gmat 0.2.0

In this version the functionality of `gmat` has been significantly extended.
In addition to the functionality already available for sampling covariance
matrices, possibly with a zero pattern specified by an undirected graph, now the
package also allows sampling correlation matrices with zero entries on their
Cholesky factor, represented by an acyclic digraph. 

## Breaking changes
Arguments for `port()` and `diagdom()` have been refactored in order to unify the
approaches for undirected graphs and the new functions for acyclic digraphs.
This has had some consequences in terms of the behaviour of the two functions.

* The arguments have been reordered taking into account their relative
  importance. In particular, `ug` is now the fourth argument, instead of the
  second one, for both functions.
* Now the function does not generate an undirected graph if parameter `ug` is
  not provided, unless explicitly stating a parameter `d < 1`.
* The default values for `p` and `d` have been changed. Now by default
  one `3 x 3` full matrix is returned.
* The argument `rentries` has been removed for both functions. In the future maybe this
  argument is reintroduced with a more complete checking of its validity
  depending on the properties of the function.
* The argument `k` has been removed from `diagdom()`, and its functionality is now
  implemented by the utility function `set_cond_number()`.

## New features
The main addition in this version are functions for correlation matrix sampling,
possibly with constraints on the upper Cholesky factorization, which correspond
to an acyclic digraph representation. Some side utility functions are also
provided.

* Added functions `chol_mh()`, `chol_iid()` and `chol_polar()`. These three
  functions return a sample of correlation matrices, possibly with an average
  percentage of zeros in their upper Cholesky factor, which can be also
  predefined by a given acyclic digraph. See more details at their
  documentation.
* Functions `mh_sphere()` and `mh_u()`, used by `chol_mh()`,
  allows to sample the upper Cholesky factor of a correlation
  matrix by sampling vectors on hemispheres of different
  dimensions. More information on their documentation.
* Added utility function `rgraph()`, which is a simple wrapper of some
  functionality in package [igraph](https://CRAN.R-project.org/package=igraph)
  for random graph generation.
* New function `anti_t()` computes the anti transpose of a matrix. This is
  mainly useful for testing, since it is involved in the acyclic digraph
  representation of the upper and lower Cholesky factors.
* Utility function `vectorize()` for extracting the upper/lower triangle in a
  covariance/correlation matrix sample as returned by the functions in the
  package.

## Minor improvements
* Updated documentation and examples for `port()` and `diagdom()`.
* Now [igraph](https://CRAN.R-project.org/package=igraph) package is not
  imported into the `NAMESPACE`, but instead explicitly called throughout the
  package using `::`.
* Removed a seemingly unnecessary package registration.

