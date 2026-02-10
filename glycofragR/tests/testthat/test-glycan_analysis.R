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

test_that("glycan_analysis_theoretical_df returns native R dataframe", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N", fragment_types = c("BY"), max_structures = 10)
  df <- glycan_analysis_theoretical_df(analysis)
  expect_s3_class(df, "data.frame")
  # Verify it's native R, not Python object
  expect_false(inherits(df, "python.builtin.dict"))
})

test_that("glycan_analysis generates fragments with BY series correctly", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  # Test with glycan code that should generate many fragments
  analysis <- glycan_analysis_new("6632", "N", 
                                   fragment_types = c("BY"),
                                   max_structures = 5,
                                   charges = c(1, 2))
  df <- glycan_analysis_theoretical_df(analysis)
  
  # Should generate more than just 3 diagnostic ions
  expect_true(nrow(df) > 10, 
              label = sprintf("Expected >10 fragments, got %d", nrow(df)))
})

test_that("glycan_analysis column names are preserved", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N", 
                                   fragment_types = c("BY"),
                                   max_structures = 5)
  df <- glycan_analysis_theoretical_df(analysis)
  
  # Check for mass-related columns
  mass_cols <- names(df)[grepl("[Mm]ass", names(df))]
  expect_true(length(mass_cols) > 0,
              label = sprintf("Mass column not found. Columns: %s", paste(names(df), collapse=", ")))
})

test_that("glycan_analysis_summary returns data.frame", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N")
  summary <- glycan_analysis_summary(analysis)
  expect_s3_class(summary, "data.frame")
  expect_true(nrow(summary) > 0)
})

test_that("glycan_analysis_clean_df returns deduplicated fragments", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  analysis <- glycan_analysis_new("4501", "N",
                                   fragment_types = c("BY"),
                                   max_structures = 10)
  
  theoretical <- glycan_analysis_theoretical_df(analysis)
  clean <- glycan_analysis_clean_df(analysis)
  
  # Clean should have same or fewer rows (deduplicated)
  expect_true(nrow(clean) <= nrow(theoretical),
              label = sprintf("Clean has %d rows, theoretical has %d", nrow(clean), nrow(theoretical)))
})
