# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes \
#     --errors-for-leak-kinds=definite --suppressions=suppressions.supp \
#     --error-exitcode=1" --vanilla < memcheck/tests.R
# First build and install the package.
testthat::test_local()
