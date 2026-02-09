test_that("peptide_new creates object", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  peptide <- peptide_new("PEPTIDE")
  expect_true(!is.null(peptide))
})

test_that("peptide_fragments returns fragments", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  peptide <- peptide_new("PEPTIDE")
  frags <- peptide_fragments(peptide, charge_states = c(1, 2))
  expect_type(frags, "list")
  expect_true(length(frags) > 0)
})

test_that("peptide with modifications works", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  peptide <- peptide_new("LCPDCPLLAPLNDSR", mod_string = "C:CAM")
  frags <- peptide_fragments(peptide)
  expect_type(frags, "list")
})
