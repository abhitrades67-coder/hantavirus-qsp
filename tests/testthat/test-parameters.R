test_that("get_parameters() returns a finite, named, scalar-numeric list", {
  p <- get_parameters()
  expect_type(p, "list")
  expect_gt(length(p), 100)
  expect_true(all(nzchar(names(p))))
  ok <- vapply(p, function(x) is.numeric(x) && length(x) == 1L && is.finite(x),
               logical(1))
  expect_true(all(ok),
              info = paste("non-scalar/non-finite:",
                           paste(names(p)[!ok], collapse = ", ")))
})

test_that("build_parameter_table() values match get_parameters() exactly", {
  # Regression guard: the table's value column must be sourced from
  # get_parameters(), never from the legacy value= literals in the row builders.
  p   <- get_parameters()
  tbl <- build_parameter_table(p)
  shared <- intersect(tbl$parameter, names(p))
  expect_gt(length(shared), 100)
  tbl_vals <- tbl$value[match(shared, tbl$parameter)]
  par_vals <- vapply(shared, function(nm) as.numeric(p[[nm]]), numeric(1))
  expect_equal(tbl_vals, unname(par_vals))
})

test_that("every model parameter is documented in the parameter table", {
  p   <- get_parameters()
  tbl <- build_parameter_table(p)
  expect_true(all(names(p) %in% tbl$parameter),
              info = paste("undocumented:",
                           paste(setdiff(names(p), tbl$parameter), collapse = ", ")))
})
