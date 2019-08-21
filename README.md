<!-- README.md is generated from README.Rmd. Please edit that file -->

# gmat

An R package for simulating positive definite matrices constrained by
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

The only R package required for `gmat` is `igraph`, which can also be
installed from CRAN.

## An example of use

First, we generate a random undirected graph with `3` nodes and density
`0.5`. Then we generate, using our `port()` function, `2` matrices
consistent with such random graphical structure.

``` r
library(gmat)

ug <- rgraph(p = 3, d = 0.5)
igraph::print.igraph(ug)
#> IGRAPH f4e129f U--- 3 1 -- Erdos renyi (gnp) graph
#> + attr: name (g/c), type (g/c), loops (g/l), p (g/n)
#> + edge from f4e129f:
#> [1] 2--3
port(N = 2, ug = ug)
#> , , 1
#> 
#>      [,1]      [,2]      [,3]
#> [1,]    1 0.0000000 0.0000000
#> [2,]    0 1.0000000 0.1668523
#> [3,]    0 0.1668523 1.0000000
#> 
#> , , 2
#> 
#>      [,1]      [,2]      [,3]
#> [1,]    1 0.0000000 0.0000000
#> [2,]    0 1.0000000 0.9987067
#> [3,]    0 0.9987067 1.0000000
```

We appreciate how the zero pattern is shared by all of the simulated
matrices. The return value is an array, and so the individual matrices
can be accessed as `matrices[, , n]`, where `n` is the index of the
matrix we want to retrieve from the sample, ranging from `1` to `N`.

We may also sample correlation matrices using i.i.d. coefficients in
their upper Cholesky factor `U`.

``` r
chol_iid(N = 2)
#> , , 1
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.3582745 -0.8077264
#> [2,]  0.3582745  1.0000000 -0.6780165
#> [3,] -0.8077264 -0.6780165  1.0000000
#> 
#> , , 2
#> 
#>            [,1]        [,2]        [,3]
#> [1,]  1.0000000 -0.77302566 -0.31058295
#> [2,] -0.7730257  1.00000000 -0.09241877
#> [3,] -0.3105829 -0.09241877  1.00000000
```

A specific zero pattern can be enforced in `U` using an acyclic digraph.

``` r
dag <- rgraph(p = 3, d = 0.5, dag = TRUE)
m <- chol_iid(dag = dag)[, , 1]
L <- t(chol(anti_t(m)))
U <- t(anti_t(L))
igraph::print.igraph(dag)
#> IGRAPH 167ee94 D--- 3 1 -- 
#> + edge from 167ee94:
#> [1] 1->2
print(U)
#>          [,1]       [,2] [,3]
#> [1,] 0.764341 -0.6448122    0
#> [2,] 0.000000  1.0000000    0
#> [3,] 0.000000  0.0000000    1
print(m)
#>            [,1]       [,2] [,3]
#> [1,]  1.0000000 -0.6448122    0
#> [2,] -0.6448122  1.0000000    0
#> [3,]  0.0000000  0.0000000    1
```

See more examples and paper references at [the documentation
website](https://irenecrsn.github.io/gmat/) for the package.
