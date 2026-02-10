test_that("peptide creates object", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  pep <- peptide("PEPTIDE")
  expect_true(!is.null(pep))
})

test_that("peptide_fragments returns list of fragments", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  pep <- peptide("PEPTIDE")
  frags <- peptide_fragments(pep, charge_states = c(1, 2))
  expect_type(frags, "list")
  expect_true(length(frags) > 0)
  # Should have multiple fragment types
  expect_true(length(names(frags)) > 0)
})

test_that("peptide with modifications works", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  pep <- peptide("LCPDCPLLAPLNDSR", mod_string = "C:CAM")
  frags <- peptide_fragments(pep, charge_states = c(1, 2))
  expect_type(frags, "list")
  expect_true(length(frags) > 0)
})

test_that("peptide handles sequence correctly", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  pep <- peptide("MSLSPDGK")
  expect_true(!is.null(pep))
  frags <- peptide_fragments(pep, charge_states = c(1))
  expect_type(frags, "list")
  # Should generate fragments (b and y ions)
  expect_true(length(frags) > 0)
})
