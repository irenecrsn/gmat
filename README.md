<!-- README.md is generated from README.Rmd. Please edit that file -->
gmat
====

An R package for simulating positive definite matrices constrained by acyclic directed and undirected graphs.

[![Build Status](https://travis-ci.org/irenecrsn/gmat.svg?branch=master)](https://travis-ci.org/irenecrsn/gmat) [![codecov](https://codecov.io/gh/irenecrsn/gmat/branch/dev/graph/badge.svg)](https://codecov.io/gh/irenecrsn/gmat) [![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](http://cran.r-project.org/package=gmat) [![CRAN status](http://www.r-pkg.org/badges/version/gmat)](http://cran.r-project.org/package=gmat)

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
#> IGRAPH 39c3a3f U--- 3 3 -- Erdos renyi (gnp) graph
#> + attr: name (g/c), type (g/c), loops (g/l), p (g/n)
#> + edges from 39c3a3f:
#> [1] 1--2 1--3 2--3
port(N = 2, ug = ug)
#> , , 1
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 1.7068774 0.9014150 0.6500106
#> [2,] 0.9014150 0.6370322 0.3225905
#> [3,] 0.6500106 0.3225905 0.8021840
#> 
#> , , 2
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 1.4033656 1.4240124 0.8930618
#> [2,] 1.4240124 1.6007928 0.9213161
#> [3,] 0.8930618 0.9213161 1.0317726
```

We appreciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.

We may also sample correlation matrices using i.i.d. coefficients in their upper Cholesky factor `U`.

``` r
chol_iid(N = 2)
#> , , 1
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.2145312 -0.1186931
#> [2,] -0.2145312  1.0000000 -0.1950826
#> [3,] -0.1186931 -0.1950826  1.0000000
#> 
#> , , 2
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1199003 -0.3526756
#> [2,] -0.1199003  1.0000000 -0.5379672
#> [3,] -0.3526756 -0.5379672  1.0000000
```

A specific zero pattern can be enforced in `U` using an acyclic digraph.

``` r
dag <- rgraph(p = 3, d = 0.5, dag = TRUE)
m <- chol_iid(dag = dag)[, , 1]
L <- t(chol(anti_t(m)))
U <- t(anti_t(L))
igraph::print.igraph(dag)
#> IGRAPH dde1742 D--- 3 2 -- 
#> + edges from dde1742:
#> [1] 1->3 2->3
print(U)
#>           [,1]          [,2]       [,3]
#> [1,] 0.8050665 -1.288206e-16 -0.5931845
#> [2,] 0.0000000  8.618367e-01 -0.5071859
#> [3,] 0.0000000  0.000000e+00  1.0000000
print(m)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.3008549 -0.5931845
#> [2,]  0.3008549  1.0000000 -0.5071859
#> [3,] -0.5931845 -0.5071859  1.0000000
```

See more examples and paper references at [the documentation website](https://irenecrsn.github.io/gmat/) for the package.
