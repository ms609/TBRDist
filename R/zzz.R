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
    "Have you refreshed the package meta with codemetar::write_codemeta()?",
    "Have you updated the version number in .zenodo.json, NEWS & DESCRIPTION?"
  )
}

#rhub::check_on_windows()
#rhub::check_with_rdevel() # redundifies check_on_debian()
#rhub::check_on_ubuntu()
#rhub::check_on_fedora()
#rhub::check_on_centos()
#rhub::check_with_valgrind() # runs the build and check on Linux, in valgrind to find memory leaks and pointer errors.
#rhub::check_with_sanitizers() # runs all package package tests, examples and vignettes with Address Sanitizer and Undefined Behavior Sanitizer.
#list_my_checks() # list_package_checks
