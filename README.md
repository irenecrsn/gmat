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

First, we generate a random undirected graph with `3` nodes and density `0.5`. Then we generate, using our `gmat::port` function, `2` matrices consistent with such random graphical structure.

``` r
ug <- gmat::rgraph(p = 3, d = 0.5)
igraph::print.igraph(ug)
#> IGRAPH 986892c U--- 3 2 -- Erdos renyi (gnp) graph
#> + attr: name (g/c), type (g/c), loops (g/l), p (g/n)
#> + edges from 986892c:
#> [1] 1--3 2--3
gmat::port(N = 2, ug = ug)
#> , , 1
#> 
#>          [,1]      [,2]      [,3]
#> [1,] 1.038434 0.0000000 1.1385855
#> [2,] 0.000000 0.6603003 0.5554253
#> [3,] 1.138586 0.5554253 1.8361219
#> 
#> , , 2
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 0.9460005 0.0000000 0.6489812
#> [2,] 0.0000000 0.2444352 0.0845033
#> [3,] 0.6489812 0.0845033 0.9609724
```

We appreciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.

We may also sample correlation matrices using i.i.d. coefficients in their upper Cholesky factor `U`.

``` r
gmat::chol_iid(N = 2)
#> , , 1
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1366272 -0.6421917
#> [2,] -0.1366272  1.0000000 -0.5668612
#> [3,] -0.6421917 -0.5668612  1.0000000
#> 
#> , , 2
#> 
#>            [,1]        [,2]        [,3]
#> [1,]  1.0000000 -0.45790243 -0.44317638
#> [2,] -0.4579024  1.00000000 -0.09992628
#> [3,] -0.4431764 -0.09992628  1.00000000
```

A specific zero pattern can be enforced in `U` using an acyclic digraph.

``` r
dag <- gmat::rgraph(p = 3, d = 0.5, dag = TRUE)
m <- gmat::chol_iid(dag = dag)[, , 1]
L <- t(chol(gmat::anti_t(m)))
U <- t(gmat::anti_t(L))
igraph::print.igraph(dag)
#> IGRAPH 3cbf8e0 D--- 3 1 -- 
#> + edge from 3cbf8e0:
#> [1] 1->2
print(U)
#>           [,1]      [,2] [,3]
#> [1,] 0.8979351 -0.440128    0
#> [2,] 0.0000000  1.000000    0
#> [3,] 0.0000000  0.000000    1
print(m)
#>           [,1]      [,2] [,3]
#> [1,]  1.000000 -0.440128    0
#> [2,] -0.440128  1.000000    0
#> [3,]  0.000000  0.000000    1
```

See more examples and paper references at [the documentation website](https://irenecrsn.github.io/gmat/) for the package.
