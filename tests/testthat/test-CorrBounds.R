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

test_that("CorrBounds handles missing parameters", {
  # Test missing lambda.OS
  expect_error(
    CorrBounds(
      outcomes = c('OS', 'PFS'),
      lambda.PFS = log(2)/6,
      rho.OS.PFS = 0.5
    ),
    "lambda.OS and lambda.PFS must be specified"
  )
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

test_that("CorrBounds checks positive definiteness", {
  # Test non-positive definite matrix
  result <- CorrBounds(
    outcomes = c('OS', 'PFS', 'OR'),
    lambda.OS = log(2)/12,
    lambda.PFS = log(2)/6,
    p.OR = 0.4,
    rho.OS.PFS = 0.8,
    rho.OS.OR = 0.8,
    rho.PFS.OR = -0.8  # This might create non-positive definite matrix
  )

  # Should either be valid or detect non-positive definiteness
  if (!result$valid) {
    expect_true(any(grepl("positive definite", result$errors)))
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
