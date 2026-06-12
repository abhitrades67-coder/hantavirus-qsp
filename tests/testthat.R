library(testthat)

# The Hantavirus QSP project is run by sourcing its R files directly rather than
# by installing it as a package, so the tests locate the project root and source
# the model themselves (see tests/testthat/setup.R).
#
# Run the suite from the project root with:
#   Rscript -e 'testthat::test_dir("tests/testthat")'
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
