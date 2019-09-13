# gmat 0.2.1

This update of `gmat` incorporates functions necessary to run the experiments in
[Córdoba et al. (2019)](https://arxiv.org/abs/1909.01062), and some minor
improvements and bug fixes.

## Breaking changes
* Now we also import package
  [gRbase](https://CRAN.R-project.org/package=gRbase), mainly for its
  functionality regarding triangulation and chordal graphs.
* Function `chol_polar()` has been removed because a faster implementation has
  been recently provided by the
  [randcorr](https://CRAN.R-project.org/package=randcorr) package.
* Function `anti_t()` has been replaced by a more useful function `uchol()`,
  which by transposing with respect to the antidiagonal obtains the upper
  Cholesky factor of a positive definite diagonal matrix (see Córdoba et al.
  2019, and the documentation of the function).
* All sampling functions now return correlation matrices instead of covariance
  matrix. This is more coherent with the return value of functions already added
  in version 0.2.0, and also has yielded more numerical stability and less
  execution time.

## New features
* New function `port_chol()` which implements a combination of the algorithms in
  `chol_mh()` and `port()`, see Córdoba et al. (2019) for the details.
* Utility function `ug_to_dag()` for obtaining an acyclic digraph with no
  v-structures from the chordal cover of an undirected graph.

## Minor improvements and bug fixes
* Styled R files using [styler](https://CRAN.R-project.org/package=styler)
  package.
* The user can now specify how to generate the random entries in the initial
  factors for both `port()` and `diagdom()`.
* Improved sample initialization in `port()`, now it is faster.
* Function `rgraph()` allows now to generate acyclic digraphs with any
  topological sorting of their nodes, by setting to `FALSE` its new `ordered`
  argument. Examples have also been added to the documentation of this function.
* Now function `mh_u()` works for any topological sort of a dag, not just the
  standard `1, ..., p`.

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

