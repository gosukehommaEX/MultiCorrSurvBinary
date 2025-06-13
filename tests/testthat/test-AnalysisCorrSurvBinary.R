test_that("AnalysisCorrSurvBinary performs event-driven analysis", {
  # Generate test data
  arm_params <- list(
    arm1 = list(
      mst.OS = 18, mst.PFS = 12, p.OR = 0.6, n = 100,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    ),
    arm2 = list(
      mst.OS = 12, mst.PFS = 8, p.OR = 0.4, n = 100,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    )
  )

  set.seed(123)
  test_data <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 10,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = arm_params,
    tau = 24,
    validate.bounds = FALSE
  )

  # Perform analysis
  analysis_results <- AnalysisCorrSurvBinary(
    data = test_data,
    E = c(25, 50),
    prioritize = "OS",
    subgroup.prioritize = c("entire"),
    alternative = "greater"
  )

  expect_s3_class(analysis_results, "tbl_df")
  expect_true(nrow(analysis_results) > 0)

  # Check required columns
  expected_cols <- c("sim", "analysis_event", "timing_subgroup", "analysis_subgroup",
                     "outcome", "comparison", "treatment_arm", "control_arm",
                     "analysis_time", "n_treatment", "n_control",
                     "events_treatment", "events_control", "pvalue")
  expect_true(all(expected_cols %in% names(analysis_results)))

  # Check analysis events
  expect_equal(sort(unique(analysis_results$analysis_event)), c(25, 50))

  # Check outcomes
  expect_equal(sort(unique(analysis_results$outcome)), c("OR", "OS", "PFS"))
})

test_that("AnalysisCorrSurvBinary validates inputs", {
  # Test invalid prioritize
  expect_error(
    AnalysisCorrSurvBinary(
      data = data.frame(sim = 1, ARM = "arm1", OS = 1),
      E = c(10),
      prioritize = "INVALID",
      subgroup.prioritize = c("entire"),
      alternative = "greater"
    ),
    "prioritize must be either 'OS' or 'PFS'"
  )

  # Test invalid alternative
  expect_error(
    AnalysisCorrSurvBinary(
      data = data.frame(sim = 1, ARM = "arm1", OS = 1),
      E = c(10),
      prioritize = "OS",
      subgroup.prioritize = c("entire"),
      alternative = "invalid"
    ),
    "alternative must be 'greater', 'less', or 'two.sided'"
  )
})

test_that("AnalysisCorrSurvBinary handles subgroups", {
  # Generate test data with subgroups
  arm_params_subgroups <- list(
    arm1 = list(
      sub1 = list(
        mst.OS = 18, mst.PFS = 10, p.OR = 0.5, n = 50,
        rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
      ),
      sub2 = list(
        mst.OS = 14, mst.PFS = 7, p.OR = 0.4, n = 50,
        rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
      )
    ),
    arm2 = list(
      sub1 = list(
        mst.OS = 12, mst.PFS = 5, p.OR = 0.3, n = 50,
        rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
      ),
      sub2 = list(
        mst.OS = 10, mst.PFS = 5, p.OR = 0.2, n = 50,
        rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
      )
    )
  )

  set.seed(456)
  test_data <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 5,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = arm_params_subgroups,
    tau = 18,
    validate.bounds = FALSE
  )

  # Perform subgroup analysis
  analysis_results <- AnalysisCorrSurvBinary(
    data = test_data,
    E = c(20),
    prioritize = "OS",
    subgroup.prioritize = c("entire", "sub1"),
    alternative = "greater"
  )

  expect_true(nrow(analysis_results) > 0)
  expect_equal(sort(unique(analysis_results$analysis_subgroup)), c("entire", "sub1"))
})

test_that("AnalysisCorrSurvBinary handles insufficient data", {
  # Test with data that has insufficient events
  minimal_data <- data.frame(
    sim = rep(1, 10),
    ARM = rep(c("arm1", "arm2"), each = 5),
    SUBGROUP = NA,
    patientID = 1:10,
    OS = rep(1, 10),
    PFS = rep(1, 10),
    OR = rep(0, 10),
    Accrual = rep(0, 10)
  )

  # This should run without error but may produce empty results
  result <- AnalysisCorrSurvBinary(
    data = minimal_data,
    E = c(20),  # More events than available
    prioritize = "OS",
    subgroup.prioritize = c("entire"),
    alternative = "greater"
  )

  expect_s3_class(result, "tbl_df")
  # May be empty due to insufficient events, which is acceptable
})

test_that("AnalysisCorrSurvBinary produces valid p-values", {
  # Generate simple test data
  arm_params <- list(
    arm1 = list(
      mst.OS = 15, mst.PFS = 10, p.OR = 0.5, n = 80,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    ),
    arm2 = list(
      mst.OS = 12, mst.PFS = 8, p.OR = 0.4, n = 80,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    )
  )

  set.seed(789)
  test_data <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 5,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = arm_params,
    tau = 24,
    validate.bounds = FALSE
  )

  analysis_results <- AnalysisCorrSurvBinary(
    data = test_data,
    E = c(30),
    prioritize = "OS",
    subgroup.prioritize = c("entire"),
    alternative = "greater"
  )

  if (nrow(analysis_results) > 0) {
    # P-values should be between 0 and 1
    expect_true(all(analysis_results$pvalue >= 0 & analysis_results$pvalue <= 1, na.rm = TRUE))
  }
})
