<!-- README.md is generated from README.Rmd. Please edit that file -->
gmat
====

An R package for simulating positive definite matrices constrained by acyclic directed and undirected graphs.

[![Build Status](https://travis-ci.org/irenecrsn/gmat.svg?branch=master)](https://travis-ci.org/irenecrsn/gmat) [![codecov](https://codecov.io/gh/irenecrsn/gmat/branch/master/graph/badge.svg)](https://codecov.io/gh/irenecrsn/gmat) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](http://cran.r-project.org/package=gmat) [![CRAN status](http://www.r-pkg.org/badges/version/gmat)](http://cran.r-project.org/package=gmat)

Installation
------------

The package is available on `CRAN`, to get the latest stable version use:

``` r
install.packages("gmat")
```

The current version on `CRAN` only includes the methods for undirected graph models.

Alternatively, using the R package `devtools` one may install the development version:

``` r
# install.packages("devtools")
devtools::install_github("irenecrsn/gmat")
```

The only R package required for `gmat` is `igraph`, which can also be installed from CRAN.

An example of use
-----------------

First, we generate a random graph with `3` nodes and density `0.25`, using the `igraph` package. Then we generate, using our `gmat::port` function, `3` matrices consistent with such random graphical structure.

``` r
ug <- igraph::sample_gnp(n = 3, p = 0.25)
matrices <- gmat::port(N = 3, ug = ug)
matrices
#> , , 1
#> 
#>         [,1]      [,2]      [,3]
#> [1,] 1.53998 0.0000000 0.0000000
#> [2,] 0.00000 0.3770087 0.0000000
#> [3,] 0.00000 0.0000000 0.4616244
#> 
#> , , 2
#> 
#>          [,1]      [,2]      [,3]
#> [1,] 1.313659 0.0000000 0.0000000
#> [2,] 0.000000 0.0145703 0.0000000
#> [3,] 0.000000 0.0000000 0.4890558
#> 
#> , , 3
#> 
#>          [,1]      [,2]      [,3]
#> [1,] 1.048078 0.0000000 0.0000000
#> [2,] 0.000000 0.0866317 0.0000000
#> [3,] 0.000000 0.0000000 0.1614992
```

We apprieciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.
