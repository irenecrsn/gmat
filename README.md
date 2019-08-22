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
#> IGRAPH 9d09b82 U--- 3 2 -- Erdos renyi (gnp) graph
#> + attr: name (g/c), type (g/c), loops (g/l), p (g/n)
#> + edges from 9d09b82:
#> [1] 1--3 2--3
port(N = 2, ug = ug)
#> , , 1
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  1.0000000 0.0000000 -0.7369786
#> [2,]  0.0000000 1.0000000  0.5321285
#> [3,] -0.7369786 0.5321285  1.0000000
#> 
#> , , 2
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.0000000 0.2265656
#> [2,] 0.0000000 1.0000000 0.3709148
#> [3,] 0.2265656 0.3709148 1.0000000
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
#> [1,]  1.0000000 -0.5646855 -0.3028382
#> [2,] -0.5646855  1.0000000 -0.3676984
#> [3,] -0.3028382 -0.3676984  1.0000000
#> 
#> , , 2
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.4525598 -0.4356729
#> [2,] -0.4525598  1.0000000 -0.4428669
#> [3,] -0.4356729 -0.4428669  1.0000000
```

A specific zero pattern can be enforced in `U` using an acyclic digraph.

``` r
dag <- rgraph(p = 3, d = 0.5, dag = TRUE)
igraph::print.igraph(dag)
#> IGRAPH 5e6add2 D--- 3 1 -- 
#> + edge from 5e6add2:
#> [1] 1->3
peo_sort <- rev(igraph::topo_sort(dag))
m <- chol_iid(dag = dag)[, , 1]
# We sort the matrix according to a Perfect Elimination Ordering (PEO)
m_sorted <- m[peo_sort, peo_sort]
# Transpose of Cholesky factor (lower Cholesky decomposition)
U_chol <- chol(m_sorted) 
# Upper Cholesky factor (upper Cholesky decomposition)
U <- t(U_chol[3:1, 3:1])
print(U)
#>           [,1] [,2]       [,3]
#> [1,] 0.7161313    0 -0.6979656
#> [2,] 0.0000000    1  0.0000000
#> [3,] 0.0000000    0  1.0000000
print(m)
#>            [,1] [,2]       [,3]
#> [1,]  1.0000000    0 -0.6979656
#> [2,]  0.0000000    1  0.0000000
#> [3,] -0.6979656    0  1.0000000
```

See more examples and paper references at [the documentation
website](https://irenecrsn.github.io/gmat/) for the package.
