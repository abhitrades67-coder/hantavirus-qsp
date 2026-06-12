# Locate the project root (the directory containing both DESCRIPTION and model/),
# regardless of whether the suite is launched from the project root or from
# tests/testthat/, then source the model so the tests can call it directly.
.find_qsp_root <- function() {
  starts <- unique(c(
    getwd(),
    normalizePath(".",      winslash = "/", mustWork = FALSE),
    normalizePath("../..",  winslash = "/", mustWork = FALSE)
  ))
  for (start in starts) {
    d <- start
    for (i in seq_len(8)) {
      if (file.exists(file.path(d, "DESCRIPTION")) &&
          dir.exists(file.path(d, "model"))) {
        return(d)
      }
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }
  stop("Could not locate the Hantavirus QSP project root (DESCRIPTION + model/).")
}

qsp_root <- .find_qsp_root()

source(file.path(qsp_root, "model", "parameters.R"))
source(file.path(qsp_root, "model", "hantavirus_qsp.R"))
source(file.path(qsp_root, "R", "pd_models.R"))
