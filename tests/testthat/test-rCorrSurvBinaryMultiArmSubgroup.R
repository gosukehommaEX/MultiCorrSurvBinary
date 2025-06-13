test_that("rCorrSurvBinaryMultiArmSubgroup handles multiple arms", {
  arm_params <- list(
    arm1 = list(
      mst.OS = 18, mst.PFS = 12, p.OR = 0.6, n = 50,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    ),
    arm2 = list(
      mst.OS = 12, mst.PFS = 8, p.OR = 0.4, n = 50,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    )
  )

  set.seed(123)
  result <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 5,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = arm_params,
    tau = 24,
    validate.bounds = FALSE
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 5 * 100)  # nsim * total n
  expect_named(result, c("sim", "ARM", "SUBGROUP", "patientID", "OS", "PFS", "OR", "Accrual"))
  expect_equal(sort(unique(result$ARM)), c("arm1", "arm2"))
  expect_true(all(is.na(result$SUBGROUP)))  # No subgroups in this example
})

test_that("rCorrSurvBinaryMultiArmSubgroup handles subgroups", {
  arm_params_subgroups <- list(
    arm1 = list(
      sub1 = list(
        mst.OS = 18, mst.PFS = 10, p.OR = 0.5, n = 30,
        rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
      ),
      sub2 = list(
        mst.OS = 14, mst.PFS = 7, p.OR = 0.4, n = 20,
        rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
      )
    )
  )

  set.seed(456)
  result <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 3,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = arm_params_subgroups,
    tau = 18,
    validate.bounds = FALSE
  )

  expect_equal(sort(unique(result$SUBGROUP)), c("sub1", "sub2"))
  expect_equal(nrow(result), 3 * 50)  # nsim * total n

  # Check subgroup sizes
  subgroup_sizes <- result %>%
    filter(sim == 1) %>%
    group_by(SUBGROUP) %>%
    summarise(n = n(), .groups = 'drop')

  expect_equal(subgroup_sizes$n[subgroup_sizes$SUBGROUP == "sub1"], 30)
  expect_equal(subgroup_sizes$n[subgroup_sizes$SUBGROUP == "sub2"], 20)
})

test_that("rCorrSurvBinaryMultiArmSubgroup validates parameters", {
  # Test missing required parameters
  invalid_params <- list(
    arm1 = list(
      mst.OS = 18,
      # Missing other required parameters
      n = 50
    )
  )

  expect_error(
    rCorrSurvBinaryMultiArmSubgroup(
      nsim = 1,
      outcomes = c('OS', 'PFS', 'OR'),
      arm.params = invalid_params,
      tau = 12
    ),
    "Missing required parameters"
  )
})

test_that("rCorrSurvBinaryMultiArmSubgroup validates arm.params structure", {
  expect_error(
    rCorrSurvBinaryMultiArmSubgroup(
      nsim = 1,
      outcomes = c('OS'),
      arm.params = "invalid",  # Should be a list
      tau = 12
    ),
    "arm.params must be a named list"
  )
})

test_that("rCorrSurvBinaryMultiArmSubgroup generates unique patient IDs", {
  arm_params <- list(
    arm1 = list(
      mst.OS = 18, mst.PFS = 12, p.OR = 0.6, n = 25,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    ),
    arm2 = list(
      mst.OS = 12, mst.PFS = 8, p.OR = 0.4, n = 25,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    )
  )

  set.seed(789)
  result <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 2,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = arm_params,
    tau = 24,
    validate.bounds = FALSE
  )

  # Check patient IDs are unique within each arm and simulation
  id_check <- result %>%
    group_by(sim, ARM) %>%
    summarise(
      max_id = max(patientID),
      n_unique = n_distinct(patientID),
      n_total = n(),
      .groups = 'drop'
    )

  expect_true(all(id_check$n_unique == id_check$n_total))
  expect_true(all(id_check$max_id == 25))  # Should be 1 to n for each arm
})

test_that("rCorrSurvBinaryMultiArmSubgroup handles mixed arm structures", {
  # Mix of single population and subgroups
  mixed_params <- list(
    arm1 = list(
      # Single population
      mst.OS = 18, mst.PFS = 12, p.OR = 0.6, n = 50,
      rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
    ),
    arm2 = list(
      # Subgroups
      sub1 = list(
        mst.OS = 12, mst.PFS = 8, p.OR = 0.4, n = 30,
        rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
      ),
      sub2 = list(
        mst.OS = 10, mst.PFS = 6, p.OR = 0.3, n = 20,
        rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
      )
    )
  )

  set.seed(101112)
  result <- rCorrSurvBinaryMultiArmSubgroup(
    nsim = 2,
    outcomes = c('OS', 'PFS', 'OR'),
    arm.params = mixed_params,
    tau = 24,
    validate.bounds = FALSE
  )

  # Check structure
  arm1_data <- result %>% filter(ARM == "arm1")
  arm2_data <- result %>% filter(ARM == "arm2")

  expect_true(all(is.na(arm1_data$SUBGROUP)))  # No subgroups for arm1
  expect_false(all(is.na(arm2_data$SUBGROUP)))  # Has subgroups for arm2
  expect_equal(sort(unique(arm2_data$SUBGROUP)), c("sub1", "sub2"))
})
