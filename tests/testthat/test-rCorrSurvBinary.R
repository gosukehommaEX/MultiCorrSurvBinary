test_that("rCorrSurvBinary generates correct data structure", {
  set.seed(123)
  result <- rCorrSurvBinary(
    nsim = 5,
    outcomes = c('OS', 'PFS', 'OR'),
    n = 50,
    mst.OS = 12,
    mst.PFS = 6,
    p.OR = 0.4,
    rho.OS.PFS = 0.5,
    rho.OS.OR = 0.3,
    rho.PFS.OR = 0.4,
    tau = 24,
    seed = 123
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 5 * 50)  # nsim * n
  # Check that required columns exist (order may vary)
  expected_cols <- c("sim", "patientID", "OS", "PFS", "OR", "Accrual")
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(unique(result$sim), 1:5)
  expect_true(all(result$patientID >= 1 & result$patientID <= 50))
})

test_that("rCorrSurvBinary enforces OS >= PFS constraint", {
  set.seed(456)
  result <- rCorrSurvBinary(
    nsim = 10,
    outcomes = c('OS', 'PFS'),
    n = 100,
    mst.OS = 12,
    mst.PFS = 6,
    rho.OS.PFS = 0.7,
    tau = 24
  )

  # Check constraint is satisfied
  expect_true(all(result$OS >= result$PFS))
})

test_that("rCorrSurvBinary validates required parameters", {
  expect_error(
    rCorrSurvBinary(
      nsim = 1,
      outcomes = c('OS'),
      n = 10,
      # Missing mst.OS
      tau = 12
    ),
    "mst.OS must be specified"
  )

  expect_error(
    rCorrSurvBinary(
      nsim = 1,
      outcomes = c('OR'),
      n = 10,
      # Missing p.OR
      tau = 12
    ),
    "p.OR must be specified"
  )
})

test_that("rCorrSurvBinary handles single outcome", {
  set.seed(789)
  result <- rCorrSurvBinary(
    nsim = 3,
    outcomes = c('OS'),
    n = 20,
    mst.OS = 15,
    tau = 12
  )

  expected_cols <- c("sim", "patientID", "OS", "Accrual")
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(nrow(result), 3 * 20)
})

test_that("rCorrSurvBinary preserves marginal distributions approximately", {
  set.seed(101112)
  result <- rCorrSurvBinary(
    nsim = 200,  # Increase sample size for better approximation
    outcomes = c('OS', 'PFS', 'OR'),
    n = 200,
    mst.OS = 12,
    mst.PFS = 6,
    p.OR = 0.4,
    rho.OS.PFS = 0.5,
    rho.OS.OR = 0.3,
    rho.PFS.OR = 0.4,
    tau = 24,
    validate.bounds = FALSE  # Skip validation for speed
  )

  # Check OS median (should be close to 12, allow for statistical variation)
  os_median <- median(result$OS)
  expect_gt(os_median, 9)    # Allow wider tolerance
  expect_lt(os_median, 15)

  # Check PFS median (should be close to 6, allow for statistical variation)
  pfs_median <- median(result$PFS)
  expect_gt(pfs_median, 4.5)  # Allow wider tolerance
  expect_lt(pfs_median, 7.5)

  # Check OR probability (should be close to 0.4)
  or_rate <- mean(result$OR)
  expect_gt(or_rate, 0.35)
  expect_lt(or_rate, 0.45)
})

test_that("rCorrSurvBinary respects prioritize parameter", {
  set.seed(131415)

  # Test prioritize = "OS"
  result_os <- rCorrSurvBinary(
    nsim = 10,
    outcomes = c('OS', 'PFS'),
    n = 100,
    mst.OS = 12,
    mst.PFS = 6,
    rho.OS.PFS = 0.8,
    prioritize = "OS",
    validate.bounds = FALSE
  )

  # Test prioritize = "PFS"
  result_pfs <- rCorrSurvBinary(
    nsim = 10,
    outcomes = c('OS', 'PFS'),
    n = 100,
    mst.OS = 12,
    mst.PFS = 6,
    rho.OS.PFS = 0.8,
    prioritize = "PFS",
    validate.bounds = FALSE
  )

  # Both should satisfy OS >= PFS constraint
  expect_true(all(result_os$OS >= result_os$PFS))
  expect_true(all(result_pfs$OS >= result_pfs$PFS))
})

test_that("rCorrSurvBinary validates outcomes parameter", {
  expect_error(
    rCorrSurvBinary(
      nsim = 1,
      outcomes = c('INVALID'),
      n = 10,
      tau = 12
    ),
    "outcomes must contain only 'OS', 'PFS', and/or 'OR'"
  )
})

test_that("rCorrSurvBinary handles correlation validation", {
  # Test with bounds validation enabled (default)
  set.seed(161718)
  result <- rCorrSurvBinary(
    nsim = 2,
    outcomes = c('OS', 'PFS'),
    n = 50,
    mst.OS = 12,
    mst.PFS = 6,
    rho.OS.PFS = 0.3,  # Conservative correlation
    tau = 12,
    validate.bounds = TRUE
  )

  expect_s3_class(result, "tbl_df")
  expect_true(all(result$OS >= result$PFS))
})
