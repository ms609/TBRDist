# uspr

[![Build Status](https://travis-ci.org/ms609/uspr.svg?branch=master)](https://travis-ci.org/ms609/uspr)
[![codecov](https://codecov.io/gh/ms609/uspr/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/uspr)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/uspr)](https://cran.r-project.org/package=uspr)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/uspr)](https://cran.r-project.org/package=uspr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3522726.svg)](http://doi.org/10.5281/zenodo.3522725)<!--[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
-->
[![Project Status: Active – – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

This package provides functions for creating, modifying and analysing 
phylogenetic trees.  It complements packages such as 
[ape](https://cran.r-project.org/package=ape),
[phangorn](https://cran.r-project.org/package=phangorn) and
[phytools](https://cran.r-project.org/package=phytools),
aiming for efficient and robust implementations of functions, typically
applied to unweighted trees (i.e. those without edge lengths).

Version 0.1.0 is a pre-release; some functionality may change in version 1.0.0,
due early 2020.

# Installation

Install and load the library from CRAN as follows:
```
install.packages('uspr')
library('uspr')
```

If you're feeling brave, you can install the development version thus:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/uspr')
```

Please note that the 'uspr' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
