# TBRDist

[![Build Status](https://travis-ci.org/ms609/TBRDist.svg?branch=master)](https://travis-ci.org/ms609/TBRDist)
[![codecov](https://codecov.io/gh/ms609/TBRDist/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TBRDist)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TBRDist)](https://cran.r-project.org/package=TBRDist)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TBRDist)](https://cran.r-project.org/package=TBRDist)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3548333.svg)](http://doi.org/10.5281/zenodo.3548333)
[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)

'TBRDist' allows the calculation of SPR, TBR and Replug distances between
unrooted phylogenetic trees from within R, using 
[algorithms](https://github.com/cwhidden/uspr) written in C++
by [Chris Whidden](https://web.cs.dal.ca/~whidden/). Whidden and Matsen (2017)
provide more information on the motivation behind this project, 
the algorithms used and their expected performance.

The uSPR distance is a natural distance metric with respect to phylogenetic tree search, as common tree search and sampling software mainly use SPR operations (or NNI operations, a subset of SPR operations). The uSPR distance is also a lower bound on the number of lateral gene transfer events required to explain the difference between a reference/species tree and a gene tree.

Note that the uSPR distance is NP-hard to compute (as are the TBR and replug
distances) and the running time of most of the algorithms used in this software
scales exponentially with the distance computed.
This version of TBRDist is aimed at trees with up to 50 leaves and uSPR
distances up to 14.

If you are interested in comparing rooted trees in terms of SPR operations,
you should use [rspr](https://github.com/cwhidden/rspr) instead. rspr is also
much more efficient and can easily handle pairs of binary rooted trees with 
200+ leaves and distances > 50.
rspr is not yet incorporated in this R package; please 
[contact the maintainer](https://github.com/ms609/TBRDist/issues/2)
if this would be useful to you.


# Installation

Install and load the library from CRAN (coming soon!) as follows:
```
install.packages('TBRDist')
library('TBRDist')
```

You can install the development version thus:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TBRDist')
```

Advanced users wishing to access the source code should note that the project
contains [submodules](https://github.blog/2016-02-01-working-with-submodules/).
To download the contents of these subdirectories, clone the project using
`git clone --recursive https://github.com/ms609/TBRDist`,
or if you've already cloned the project, run
`git submodule update --init --recursive`.  

# Usage

If you use 'uspr' in your research, please cite:

Chris Whidden and Frederick A. Matsen IV. Calculating the Unrooted Subtree-Prune-and-Regraft Distance. eprint arXiv:1511.07529. http://arxiv.org/abs/1511.07529


# Notes

Some text in this README file originates from the ['uspr' readme](https://github.com/cwhidden/uspr/blob/master/README.md).

Please note that the 'TBRDist' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
