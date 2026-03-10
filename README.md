# TBRDist

[![R-CMD-check](https://github.com/ms609/TBRDist/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ms609/TBRDist/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/ms609/TBRDist/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TBRDist)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/TBRDist)](https://cran.r-project.org/package=TBRDist)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/TBRDist)](https://cran.r-project.org/package=TBRDist)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3548332.svg)](https://doi.org/10.5281/zenodo.3548332)
[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

'TBRDist' is an R package that calculates rearrangement distances between
phylogenetic trees using C++ algorithms
by [Chris Whidden](https://web.cs.dal.ca/~whidden/):

- **Unrooted trees** (`TBRDist()`, `USPRDist()`, `ReplugDist()`): TBR, SPR,
  and Replug distances based on the algorithms of
  [Whidden and Matsen (2017)](https://arxiv.org/abs/1511.07529).
  The C++ implementation has diverged from the upstream
  [uspr](https://github.com/cwhidden/uspr) library with optimizations
  including integer-encoded tree identity (replacing Newick string
  round-tripping), precomputed SPR lookup tables for small trees, and
  reduced heap allocation in the A\* search and branch-and-bound routines,
  yielding roughly 2.5× faster computation.
- **Rooted trees** (`RSPRDist()`): Rooted SPR distances via
  [rspr](https://github.com/cwhidden/rspr), implementing the exact
  fixed-parameter algorithm (with cluster decomposition) of Whidden, Beiko,
  and Zeh (2013). Supports exact distances, a linear-time 3-approximation, and
  maximum agreement forest output.

The uSPR distance is a natural distance metric with respect to phylogenetic
tree search, as common tree search and sampling software mainly use SPR
operations (or NNI operations, a subset of SPR operations). The uSPR distance
is also a lower bound on the number of lateral gene transfer events required to
explain the difference between a reference/species tree and a gene tree.

Note that the uSPR distance is NP-hard to compute (as are the TBR and replug
distances) and the running time of most of the algorithms used in this software
scales exponentially with the distance computed.
The unrooted functions are aimed at trees with up to 50 leaves and uSPR
distances up to 14. The rooted `RSPRDist()` can handle much larger trees
(200+ leaves, distances > 50).


# Installation

Install and load the library from CRAN as follows:
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

If you use unrooted distances (`TBRDist`, `USPRDist`, `ReplugDist`) in your
research, please cite:

Chris Whidden and Frederick A. Matsen IV. Calculating the Unrooted Subtree-Prune-and-Regraft Distance. eprint arXiv:1511.07529. https://arxiv.org/abs/1511.07529

If you use rooted SPR distances (`RSPRDist`) in your research, please cite:

Chris Whidden, Robert G. Beiko, and Norbert Zeh. Fixed-Parameter Algorithms for Maximum Agreement Forests. _SIAM Journal on Computing_ 42(4):1431–1466, 2013. https://doi.org/10.1137/110845045


# Notes

Some text in this README file originates from the ['uspr' readme](https://github.com/cwhidden/uspr/blob/master/README.md).

Please note that the 'TBRDist' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
