# gmat - An R package for simulating covariance and concentration graph matrices

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
method by setting `method = domdiag` in the parameter list.

