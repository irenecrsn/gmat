# gmat

An R package for simulating correlation matrices possibly constrained by
acyclic directed and undirected graphs.

[![CRAN
status](http://www.r-pkg.org/badges/version/gmat)](https://CRAN.R-project.org/package=gmat)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/gmat)](https://CRAN.R-project.org/package=gmat)

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

The other R packages required for `gmat` are `igraph` and `gRbase`,
which can be installed from CRAN and Bioconductor.

## Overview

The package mostly implements methods described in the following papers:

  - Córdoba I., Varando G., Bielza C., Larrañaga P. **A fast
    Metropolis-Hastings method for generating random correlation
    matrices**. *Lecture Notes in Computer Science* (IDEAL 2018), vol
    11314, pp. 117-124, 2018.
  - Córdoba I., Varando G., Bielza C., Larrañaga P. **A partial
    orthogonalization method for simulating covariance and concentration
    graph matrices**. *Proceedings of Machine Learning Research* (PGM
    2018), vol 72, pp. 61-72, 2018.
  - Córdoba I., Varando G., Bielza C., Larrañaga P. **On generating random
    Gaussian graphical models**. *International Journal of Approximate
    Reasoning*, vol 125, pp. 240-250, 2020.

See examples of use and more at package's manual.
