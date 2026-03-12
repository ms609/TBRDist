# Code to be run with  
#   R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < memcheck/tests.R
# First build and install the package.
library("TBRDist")
devtools::test()
