## Test environments
* Local Ubuntu 16.04.6 LTS, R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.1, R devel and R 3.5.3
* Mac OS X 10.13.3 (on travis-ci), R 3.6.1 and R 3.5.3
* win-builder (R 3.6.1, R devel and R 3.5.3)

## R CMD check results
There were no ERRORs or WARNINGs.

On win-builder, in all R versions, there is a NOTE informing about the
maintainer and indicating possibly misspelled words in DESCRIPTION (Córdoba, et,
al), which are all false positives.

## Downstream dependencies
There are currently no downstream dependencies for this package.
