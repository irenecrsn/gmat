
<!-- README.md is generated from README.Rmd. Please edit that file -->
gmat
====

An R package for simulating covariance and concentration graph matrices.

[![Build Status](https://travis-ci.org/irenecrsn/gmat.svg?branch=master)](https://travis-ci.org/irenecrsn/gmat) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](http://cran.r-project.org/package=gmat) [![CRAN status](http://www.r-pkg.org/badges/version/gmat)](http://cran.r-project.org/package=gmat)

This package implements the methods described in the [paper](http://proceedings.mlr.press/v72/cordoba18a.html):

> Córdoba, I., Varando, G., Bielza, C. and Larrañaga, P. A partial orthogonalization method for simulating covariance and concentration graph matrices, Proceedings of Machine Learning Research (PGM 2018), vol. 72, pp. 61 - 72, 2018.

It is governed by two main functions, `port` and `diagdom`, which take as input an undirected graph structure `ug` and a sample size `N`, and generates `N` covariance/concentration matrices subject to the zeros imposed by the undirected graph structure `ug`.

The sample is generated by `port` with the method described in the above paper. The classical diagonal dominance method is the one used by function `diagdom`.

Installation
------------

Using the R package `devtools` one may install the development version from this branch:

``` r
# install.packages("devtools")
devtools::install_github("irenecrsn/gmat", ref = "rchol")
```

The only R package required for `gmat` is `igraph`, which can be installed from CRAN.

An example of use
-----------------

First, we generate a random graph with `3` nodes and density `0.25`, using the `igraph` package. Then we generate, using our `gmat::port` function, `3` matrices consistent with such random graphical structure.

``` r
ug <- igraph::sample_gnp(n = 3, p = 0.25)
matrices <- gmat::port(N = 3, ug = ug)
matrices
#> , , 1
#> 
#>          [,1]      [,2]      [,3]
#> [1,] 1.232688 0.0000000 0.0000000
#> [2,] 0.000000 0.1412365 0.2428490
#> [3,] 0.000000 0.2428490 0.4274721
#> 
#> , , 2
#> 
#>          [,1]       [,2]       [,3]
#> [1,] 1.172516  0.0000000  0.0000000
#> [2,] 0.000000  0.0557089 -0.0297474
#> [3,] 0.000000 -0.0297474  0.0340088
#> 
#> , , 3
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 0.9762841 0.0000000 0.0000000
#> [2,] 0.0000000 0.8823712 0.4775285
#> [3,] 0.0000000 0.4775285 0.2587886
```

We apprieciate how the zero pattern is shared by all of the simulated matrices. The return value is an array, and so the individual matrices can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we want to retrieve from the sample, ranging from `1` to `N`.
