test_that("glycopeptide_new creates object", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
  expect_true(!is.null(gp))
})

test_that("glycopeptide_fragments returns fragments", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
  frags <- glycopeptide_fragments(gp)
  expect_type(frags, "list")
  expect_true(length(frags) > 0)
})

test_that("glycopeptide_summary returns info", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
  summary <- glycopeptide_summary(gp)
  expect_type(summary, "list")
  expect_equal(summary$peptide_sequence, "LCPDCPLLAPLNDSR")
  expect_equal(summary$glycan_code, "4501")
  expect_equal(summary$glycosylation_site, 12)
})
