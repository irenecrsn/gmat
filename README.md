<!-- README.md is generated from README.Rmd. Please edit that file -->

# gmat

An R package for simulating correlation matrices possibly constrained by
acyclic directed and undirected graphs.

[![Build
Status](https://travis-ci.com/irenecrsn/gmat.svg?branch=master)](https://travis-ci.com/irenecrsn/gmat)
[![codecov](https://codecov.io/gh/irenecrsn/gmat/branch/dev/graph/badge.svg)](https://codecov.io/gh/irenecrsn/gmat)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](https://CRAN.R-project.org/package=gmat)
[![CRAN
status](http://www.r-pkg.org/badges/version/gmat)](https://CRAN.R-project.org/package=gmat)

## Installation

The package is available on `CRAN`, to get the latest stable version
use:

``` r
install.packages("gmat")
```

Alternatively, using the R package `devtools` one may install the
development version:

``` r
# install.packages("devtools")
devtools::install_github("irenecrsn/gmat")
```

The other R packages required for `gmat` are `igraph` and `gRbase`,
which can also be installed from CRAN.

## Overview

The package mostly implements methods described in the following papers:

  - Córdoba I., Varando G., Bielza C., Larrañaga P. **A fast
    Metropolis-Hastings method for generating random correlation
    matrices**. *Lecture Notes in Computer Science* (IDEAL 2018), vol
    11314, pp. 117-124, 2018.
  - Córdoba I., Varando G., Bielza C., Larrañaga P. **A partial
    orthogonalization method for simulating covariance and concentration
    graph matrices**. *Proceedings of Machine Learning Research* (PGM
    2018), vol 72, pp. 61-72,
    2018. 
  - Córdoba I., Varando G., Bielza C., Larrañaga P. **Generating random
    Gaussian graphical models**. *arXiv:1909.01062*, 2019.

## An example of use

First, we generate a random undirected graph with `3` nodes and density
`0.5`. Then we generate, using our `port()` function, `2` correlation
matrices consistent with such random graphical structure.

``` r
ug <- gmat::rgraph(p = 3, d = 0.5)
igraph::print.igraph(ug)
#> IGRAPH 564efb2 U--- 3 1 -- Erdos renyi (gnp) graph
#> + attr: name (g/c), type (g/c), loops (g/l), p (g/n)
#> + edge from 564efb2:
#> [1] 2--3
gmat::port(N = 2, ug = ug)
#> , , 1
#> 
#>      [,1]      [,2]      [,3]
#> [1,]    1 0.0000000 0.0000000
#> [2,]    0 1.0000000 0.8574111
#> [3,]    0 0.8574111 1.0000000
#> 
#> , , 2
#> 
#>      [,1]      [,2]      [,3]
#> [1,]    1  0.000000  0.000000
#> [2,]    0  1.000000 -0.642288
#> [3,]    0 -0.642288  1.000000
```

We appreciate how the zero pattern is shared by all of the simulated
matrices. The return value is an array, and so the individual matrices
can be accessed as `matrices[, , n]`, where `n` is the index of the
matrix we want to retrieve from the sample, ranging from `1` to `N`.

We may also sample correlation matrices using i.i.d. coefficients in
their upper Cholesky factor `U`.

``` r
gmat::chol_iid(N = 2)
#> , , 1
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1440407 -0.1535844
#> [2,] -0.1440407  1.0000000 -0.7167654
#> [3,] -0.1535844 -0.7167654  1.0000000
#> 
#> , , 2
#> 
#>            [,1]        [,2]        [,3]
#> [1,]  1.0000000 -0.35987932 -0.35787183
#> [2,] -0.3598793  1.00000000 -0.09487107
#> [3,] -0.3578718 -0.09487107  1.00000000
```

A specific zero pattern can be enforced in `U` using an acyclic digraph.

``` r
dag <- gmat::rgraph(p = 3, d = 0.5, dag = TRUE)
igraph::print.igraph(dag)
#> IGRAPH af3d6f6 D--- 3 1 -- 
#> + edge from af3d6f6:
#> [1] 2->3
top_sort <- igraph::topo_sort(dag)
peo_sort <- rev(top_sort)
inv_top <- order(top_sort)
inv_peo <- order(peo_sort)
m <- gmat::chol_iid(dag = dag)[, , 1]
# We sort the matrix according to both a topological and a perfect elimination
# ordering
m_top_sort <- m[top_sort, top_sort]
m_peo_sort <- m[peo_sort, peo_sort]
# Transpose of Cholesky factor (lower Cholesky decomposition)
U_chol <- chol(m_peo_sort) 
# Upper Cholesky factor (upper Cholesky decomposition)
U <- gmat::uchol(m_top_sort)
print(U_chol[inv_peo, inv_peo])
#>      [,1]       [,2] [,3]
#> [1,]    1  0.0000000    0
#> [2,]    0  0.9658269    0
#> [3,]    0 -0.2591879    1
print(U[inv_top, inv_top])
#>      [,1]      [,2]       [,3]
#> [1,]    1 0.0000000  0.0000000
#> [2,]    0 0.9658269 -0.2591879
#> [3,]    0 0.0000000  1.0000000
print(m)
#>      [,1]       [,2]       [,3]
#> [1,]    1  0.0000000  0.0000000
#> [2,]    0  1.0000000 -0.2591879
#> [3,]    0 -0.2591879  1.0000000
```

The zeros are correctly reflected in the upper Cholesky factor when
using a topological sorting, whereas they are reflected in the standard
lower Cholesky decomposition factor if we sort the matrix with a perfect
elimination ordering.

See more examples and paper references at [the documentation
website](https://irenecrsn.github.io/gmat/) for the package.
