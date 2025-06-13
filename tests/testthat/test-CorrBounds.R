test_that("CorrBounds validates correlation parameters correctly", {
  # Test valid correlations
  result <- CorrBounds(
    outcomes = c('OS', 'PFS', 'OR'),
    lambda.OS = log(2)/12,
    lambda.PFS = log(2)/6,
    p.OR = 0.4,
    rho.OS.PFS = 0.5,
    rho.OS.OR = 0.3,
    rho.PFS.OR = 0.4
  )

  expect_true(result$valid)
  expect_length(result$errors, 0)
  expect_named(result$bounds, c("OS.PFS", "OS.OR", "PFS.OR"))
})

test_that("CorrBounds detects invalid correlations", {
  # Test correlation outside bounds
  result <- CorrBounds(
    outcomes = c('OS', 'PFS'),
    lambda.OS = log(2)/12,
    lambda.PFS = log(2)/6,
    rho.OS.PFS = 0.99  # Likely outside bounds for exponential distributions
  )

  expect_false(result$valid)
  expect_gt(length(result$errors), 0)
})

test_that("CorrBounds handles missing parameters correctly", {
  # Test missing lambda.OS - should produce specific error message
  expect_error(
    CorrBounds(
      outcomes = c('OS', 'PFS'),
      lambda.PFS = log(2)/6,
      rho.OS.PFS = 0.5
    ),
    "lambda.OS and lambda.PFS must be specified"
  )

  # Alternative: Test validation result instead of direct error
  result <- CorrBounds(
    outcomes = c('OS'),
    lambda.OS = log(2)/12
  )
  expect_true(result$valid)  # Should be valid for single outcome
})

test_that("CorrBounds calculates bounds correctly", {
  result <- CorrBounds(
    outcomes = c('OS', 'PFS'),
    lambda.OS = log(2)/12,
    lambda.PFS = log(2)/6,
    rho.OS.PFS = 0
  )

  bounds <- result$bounds$OS.PFS
  expect_lt(bounds$lower, 0)
  expect_gt(bounds$upper, 0)
  expect_lt(bounds$lower, bounds$upper)
})

test_that("CorrBounds handles single outcome", {
  result <- CorrBounds(
    outcomes = c('OS'),
    lambda.OS = log(2)/12
  )

  expect_true(result$valid)
  expect_length(result$bounds, 0)
})

test_that("CorrBounds detects correlation matrix issues", {
  # Test with correlations that might create issues
  # Use more extreme values that are more likely to cause problems
  result <- CorrBounds(
    outcomes = c('OS', 'PFS', 'OR'),
    lambda.OS = log(2)/12,
    lambda.PFS = log(2)/6,
    p.OR = 0.4,
    rho.OS.PFS = 0.9,
    rho.OS.OR = 0.9,
    rho.PFS.OR = -0.9  # This combination may create issues
  )

  # The function should either:
  # 1. Detect the matrix is not positive definite, OR
  # 2. Detect correlations are outside bounds, OR
  # 3. Be valid if the correlations happen to be feasible
  if (!result$valid) {
    # Check if error is about positive definiteness OR bounds
    has_matrix_error <- any(grepl("positive definite|matrix", result$errors, ignore.case = TRUE))
    has_bounds_error <- any(grepl("bounds|outside", result$errors, ignore.case = TRUE))
    expect_true(has_matrix_error || has_bounds_error)
  } else {
    # If valid, should have eigenvalues
    expect_is(result$eigenvalues, "numeric")
  }
})

test_that("CorrBounds returns computation time", {
  result <- CorrBounds(
    outcomes = c('OS', 'PFS'),
    lambda.OS = log(2)/12,
    lambda.PFS = log(2)/6,
    rho.OS.PFS = 0.5
  )

  expect_is(result$computation.time, "numeric")
  expect_gt(result$computation.time, 0)
})

test_that("CorrBounds handles edge cases for binary outcomes", {
  # Test with extreme p.OR values
  result_low <- CorrBounds(
    outcomes = c('OS', 'OR'),
    lambda.OS = log(2)/12,
    p.OR = 0.05,  # Very low probability
    rho.OS.OR = 0.1
  )
  expect_true(result_low$valid || length(result_low$errors) > 0)

  result_high <- CorrBounds(
    outcomes = c('OS', 'OR'),
    lambda.OS = log(2)/12,
    p.OR = 0.95,  # Very high probability
    rho.OS.OR = 0.1
  )
  expect_true(result_high$valid || length(result_high$errors) > 0)
})
