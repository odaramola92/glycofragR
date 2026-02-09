test_that("set_glycofrag_conda_env works", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  expect_silent(set_glycofrag_conda_env("glycofrag"))
})

test_that("list_supported_modifications returns data", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  mods <- list_supported_modifications(verbose = FALSE)
  expect_type(mods, "list")
})
