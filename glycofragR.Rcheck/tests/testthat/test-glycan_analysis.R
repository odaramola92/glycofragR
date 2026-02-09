test_that("glycan_analysis_new creates analysis", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N")
  expect_true(!is.null(analysis))
})

test_that("glycan_analysis_structures returns structures", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N")
  structures <- glycan_analysis_structures(analysis)
  expect_type(structures, "list")
  expect_true(length(structures) > 0)
})

test_that("glycan_analysis_theoretical_df returns dataframe", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N")
  df <- glycan_analysis_theoretical_df(analysis)
  expect_s3_class(df, "data.frame")
})

test_that("glycan_analysis_summary returns summary", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N")
  summary <- glycan_analysis_summary(analysis)
  expect_s3_class(summary, "data.frame")
})
