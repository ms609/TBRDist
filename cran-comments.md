## Test environments
* Windows 10 on local machine, R 4.5.2
* Windows, ubuntu, macOS via [GitHub actions](https://github.com/ms609/TBRDist/actions/workflows/R-CMD-check.yml)
* win-builder, with `check_win_devel()`, R devel

## R CMD check results
There were no ERRORs or WARNINGs.

Two NOTEs:

* checking pragmas in C/C++ headers and code ... NOTE
  File which contains pragma(s) suppressing diagnostics: 'src/uspr/tbr.h'

  This pragma suppresses a `-Wnonnull` warning triggered by Boost's
  concept-check headers (`boost/graph/adjacency_list.hpp`), which use
  `((Model*)0)->~Model()`.  This is a known Boost issue (not a TBRDist bug)
  that fires on GCC 12+.  The suppression is tightly scoped to the two
  `#include` lines that pull in the affected Boost headers.

* checking compiled code ... NOTE

  This NOTE arises from a bug in R 4.5.x's `read_symbols_from_dll()` utility
  (`'length = 2' in coercion to 'logical(1)'`), which fails when a package DLL
  exports more than one symbol.  It is not related to TBRDist's code and does
  not appear on platforms where the utility functions correctly.

## Downstream dependencies
There are currently no downstream dependencies for this package.
