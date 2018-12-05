<!-- README.md is generated from README.Rmd. Please edit that file -->
gmat
====

An R package for simulating positive definite matrices constrained by acyclic directed and undirected graphs.

[![Build Status](https://travis-ci.org/irenecrsn/gmat.svg?branch=master)](https://travis-ci.org/irenecrsn/gmat) [![codecov](https://codecov.io/gh/irenecrsn/gmat/branch/dev/graph/badge.svg)](https://codecov.io/gh/irenecrsn/gmat) [![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](https://CRAN.R-project.org/package=gmat) [![CRAN status](http://www.r-pkg.org/badges/version/gmat)](https://CRAN.R-project.org/package=gmat)

Installation
------------

The package is available on `CRAN`, to get the latest stable version use:

``` r
install.packages("gmat")
```

Alternatively, using the R package `devtools` one may install the development version:

``` r
# install.packages("devtools")
devtools::install_github("irenecrsn/gmat")
```

The only R package required for `gmat` is `igraph`, which can also be installed from CRAN.

An example of use
-----------------

First, we generate a random undirected graph with `3` nodes and density `0.5`. Then we generate, using our `port()` function, `2` matrices consistent with such random graphical structure.

``` r
library(gmat)

ug <- rgraph(p = 3, d = 0.5)
igraph::print.igraph(ug)
#> IGRAPH 75df7f1 U--- 3 1 -- Erdos renyi (gnp) graph
#> + attr: name (g/c), type (g/c), loops (g/l), p (g/n)
#> + edge from 75df7f1:
#> [1] 1--3
port(N = 2, ug = ug)
#> , , 1
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 1.2329616 0.0000000 0.6513035
#> [2,] 0.0000000 0.1578579 0.0000000
#> [3,] 0.6513035 0.0000000 0.5332192
#> 
#> , , 2
#> 
#>          [,1]      [,2]     [,3]
#> [1,] 1.787082 0.0000000 1.296915
#> [2,] 0.000000 0.0465929 0.000000
#> [3,] 1.296915 0.0000000 1.214313
```

We appreciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.

We may also sample correlation matrices using i.i.d. coefficients in their upper Cholesky factor `U`.

``` r
chol_iid(N = 2)
#> , , 1
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  1.0000000  0.201973 -0.4340289
#> [2,]  0.2019730  1.000000 -0.6814360
#> [3,] -0.4340289 -0.681436  1.0000000
#> 
#> , , 2
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.3758306 -0.1389546
#> [2,] -0.3758306  1.0000000 -0.6183038
#> [3,] -0.1389546 -0.6183038  1.0000000
```

A specific zero pattern can be enforced in `U` using an acyclic digraph.

``` r
dag <- rgraph(p = 3, d = 0.5, dag = TRUE)
m <- chol_iid(dag = dag)[, , 1]
L <- t(chol(anti_t(m)))
U <- t(anti_t(L))
igraph::print.igraph(dag)
#> IGRAPH fdee8a0 D--- 3 1 -- 
#> + edge from fdee8a0:
#> [1] 2->3
print(U)
#>      [,1]      [,2]       [,3]
#> [1,]    1 0.0000000  0.0000000
#> [2,]    0 0.9836563 -0.1800561
#> [3,]    0 0.0000000  1.0000000
print(m)
#>      [,1]       [,2]       [,3]
#> [1,]    1  0.0000000  0.0000000
#> [2,]    0  1.0000000 -0.1800561
#> [3,]    0 -0.1800561  1.0000000
```

See more examples and paper references at [the documentation website](https://irenecrsn.github.io/gmat/) for the package.
