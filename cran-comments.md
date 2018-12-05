## Test environments
* Local Ubuntu 16.04.5 LTS, R 3.5.1
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.5.1, R devel and R 3.4.4
* Mac OS X 10.13.3 (on travis-ci), R 3.5.0 and R 3.4.4
* win-builder (R 3.5.1, R devel and R 3.4.4)

## R CMD check results
There were no ERRORs or WARNINGs.

There is one NOTE informing about the maintainer. 

There is also an indication of possible misspelled words in DESCRIPTION
(acyclic, undirected), which are all false positives.

There is an additional NOTE on win-builder with R 3.4.4, about conflicting
Author and derived Authors@R field. This seems to be because the string
`https://orcid.org` is not generated for such version. We have kept it that way
because when explicitly including ORCIDs on the Author field to solve the NOTE,
then analogous NOTEs arise on win-builder with R 3.5.0 and R devel, since these
versions do prepend `https://orcid.org`.

## Downstream dependencies
There are currently no downstream dependencies for this package.

