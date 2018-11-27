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
#>           [,1]      [,2]      [,3]
#> [1,] 0.5346747 0.0000000 0.0000000
#> [2,] 0.0000000 0.5572111 0.0000000
#> [3,] 0.0000000 0.0000000 0.5455936
#> 
#> , , 2
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 0.8960595 0.0000000 0.0000000
#> [2,] 0.0000000 0.1647816 0.0000000
#> [3,] 0.0000000 0.0000000 0.7255057
#> 
#> , , 3
#> 
#>          [,1]      [,2]      [,3]
#> [1,] 0.730354 0.0000000 0.0000000
#> [2,] 0.000000 0.3280893 0.0000000
#> [3,] 0.000000 0.0000000 0.2511559
```

We apprieciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.
