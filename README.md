# gmat - An R package for simulating covariance and concentration graph matrices
[![Build Status](https://travis-ci.org/irenecrsn/gmat.svg?branch=master)](https://travis-ci.org/irenecrsn/gmat)

## Requirements
- R packages `igraph`, `Matrix`

## Installation instructions
Probably the easiest way to install the package until it is available on CRAN is
to use the R package `devtools`:
```R
library("devtools")
install_github(repo = "irenecrsn/gmat")
```

## Main functionality
This package implements the methods described in the paper:
> Córdoba, I., Varando, G., Bielza, C. and Larrañaga, P.
> A partial orthogonalization method for simulating covariance and concentration graph matrices
> PGM 2018, Accepted.

It is governed by one main function, `rgmn`, which takes as input an undirected
graph structure `ug` and a sample size `N`, and generates `N`
covariance/concentration matrices subject to the zeros imposed by the undirected
graph structure `ug`. 

By default, it uses the method described in the above
paper; however, it is possible to also use the classical diagonal dominance
method by setting `method = "diagdom"` in the parameter list.

## An example of use

First, we generate a random graph with `3` nodes and density `0.25`, using the
`igraph` package. Then we generate, using our `gmat::rgmn` function, `3` matrices
consistent with such random graphical structure.

```R
> ug <- igraph::sample_gnp(n = 3, p = 0.25)
> matrices <- gmat::rgmn(N = 3, ug = ug, zapzeros = TRUE)
> matrices
, , 1

         [,1]      [,2]     [,3]
[1,] 1.679735 0.0000000 1.497825
[2,] 0.000000 0.2917888 0.000000
[3,] 1.497825 0.0000000 1.336177

, , 2

          [,1]      [,2]      [,3]
[1,] 0.5227167 0.0000000 0.7605577
[2,] 0.0000000 0.0348654 0.0000000
[3,] 0.7605577 0.0000000 1.1388367

, , 3

          [,1]      [,2]      [,3]
[1,] 1.8866349 0.0000000 0.9017342
[2,] 0.0000000 0.0017455 0.0000000
[3,] 0.9017342 0.0000000 0.5446746
```

The parameter `zapzeros` indicates that we want to display as zero extremely
small values. We apprieciate how the zero pattern is shared by all of the
simulated matrices. The return value is an array, and so the individual matrices
can be accessed as `matrices[, , n]`, where `n` is the index of the matrix we
want to retrieve from the sample, ranging from `1` to `N`.

