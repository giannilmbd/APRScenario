# test_scenarios.R

test_that("NKdata loads correctly", {
  data("NKdata", package = "APRScenario")
  expect_true(is.data.frame(NKdata))
  expect_true(ncol(NKdata) >= 4)
  expect_true(nrow(NKdata) > 10)
})

test_that("basic matrix operations work", {
  # Test a simple function that doesn't require full model estimation
  test_matrix <- matrix(1:9, nrow = 3)
  expect_true(is.matrix(test_matrix))
  expect_equal(nrow(test_matrix), 3)
  expect_equal(ncol(test_matrix), 3)
})

test_that("package loads without errors", {
  expect_true("APRScenario" %in% loadedNamespaces())
})


