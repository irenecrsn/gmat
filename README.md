<!-- README.md is generated from README.Rmd. Please edit that file -->
gmat
====

An R package for simulating positive definite matrices constrained by acyclic directed and undirected graphs.

[![Build Status](https://travis-ci.org/irenecrsn/gmat.svg?branch=master)](https://travis-ci.org/irenecrsn/gmat) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](http://cran.r-project.org/package=gmat) [![CRAN status](http://www.r-pkg.org/badges/version/gmat)](http://cran.r-project.org/package=gmat)

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
#> [1,] 0.9742306 0.0000000 0.2382255
#> [2,] 0.0000000 0.0174842 0.0458809
#> [3,] 0.2382255 0.0458809 0.4161402
#> 
#> , , 2
#> 
#>          [,1]       [,2]       [,3]
#> [1,] 1.315333  0.0000000  1.0398587
#> [2,] 0.000000  0.0975722 -0.0574187
#> [3,] 1.039859 -0.0574187  0.9387271
#> 
#> , , 3
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 1.5970788 0.0000000 0.8978575
#> [2,] 0.0000000 0.1401653 0.0450929
#> [3,] 0.8978575 0.0450929 0.5553658
```

We apprieciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.
