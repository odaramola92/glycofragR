test_that("glycan_new creates object", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  glycan <- glycan_new("4501", "N")
  expect_true(!is.null(glycan))
})

test_that("glycan_predict_structures returns structures", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  glycan <- glycan_new("4501", "N")
  structures <- glycan_predict_structures(glycan)
  expect_type(structures, "list")
  expect_true(length(structures) > 0)
})

test_that("glycan_generate_fragments returns fragments", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  glycan <- glycan_new("4501", "N")
  structures <- glycan_predict_structures(glycan)
  frags <- glycan_generate_fragments(glycan, structures[[1]])
  expect_type(frags, "list")
})

test_that("glycan_mass calculates mass", {
  skip_if_not(reticulate::py_available(initialize = FALSE), "Python not available")
  
  set_glycofrag_conda_env("glycofrag")
  mass <- glycan_mass("4501")
  expect_type(mass, "double")
  expect_true(mass > 0)
})
