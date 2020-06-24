#' @useDynLib TBRDist
.onUnload <- function (libpath) {
  library.dynam.unload("TBRDist", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Have you cleared GitHub issues for this release milestone?",
    "Is the code free of #TODOs?",
    "Have you checked the Vignette for sanity?",
    "Have you updated the version number in NEWS.md & DESCRIPTION?"
  )
}

# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
# check_win_devel(); rhub::check_for_cran()
# codemetar::write_codemeta()
# # revdepcheck::revdep_check()
