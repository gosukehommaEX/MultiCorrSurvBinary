test_that("TestsSurvBinary performs survival tests correctly", {
  # Create test data
  n_treatment <- 50
  n_control <- 50

  set.seed(123)
  test_data <- data.frame(
    ARM = c(rep("treatment", n_treatment), rep("control", n_control)),
    OS = c(rexp(n_treatment, rate = log(2)/15), rexp(n_control, rate = log(2)/12)),
    PFS = c(rexp(n_treatment, rate = log(2)/10), rexp(n_control, rate = log(2)/8)),
    OR = c(rbinom(n_treatment, 1, 0.6), rbinom(n_control, 1, 0.4))
  )

  # Add event indicators (assume all events occur)
  test_data$Event.OS <- 1
  test_data$Event.PFS <- 1

  # Test OS analysis
  result_os <- TestsSurvBinary(
    data = test_data,
    outcome = "OS",
    treatment_arm = "treatment",
    control_arm = "control",
    alternative = "greater"
  )

  expect_is(result_os, "list")
  expect_named(result_os, c("n_treatment", "n_control", "events_treatment", "events_control", "pvalue"))
  expect_equal(result_os$n_treatment, n_treatment)
  expect_equal(result_os$n_control, n_control)
  expect_true(result_os$pvalue >= 0 & result_os$pvalue <= 1)

  # Test PFS analysis
  result_pfs <- TestsSurvBinary(
    data = test_data,
    outcome = "PFS",
    treatment_arm = "treatment",
    control_arm = "control",
    alternative = "greater"
  )

  expect_is(result_pfs, "list")
  expect_true(result_pfs$pvalue >= 0 & result_pfs$pvalue <= 1)

  # Test OR analysis
  result_or <- TestsSurvBinary(
    data = test_data,
    outcome = "OR",
    treatment_arm = "treatment",
    control_arm = "control",
    alternative = "greater"
  )

  expect_is(result_or, "list")
  expect_true(result_or$pvalue >= 0 & result_or$pvalue <= 1)
  expect_equal(result_or$events_treatment, sum(test_data$OR[test_data$ARM == "treatment"]))
  expect_equal(result_or$events_control, sum(test_data$OR[test_data$ARM == "control"]))
})

test_that("TestsSurvBinary handles different alternative hypotheses", {
  # Create simple test data
  test_data <- data.frame(
    ARM = c(rep("treatment", 30), rep("control", 30)),
    OS = c(rexp(30, rate = log(2)/15), rexp(30, rate = log(2)/12)),
    Event.OS = 1,
    OR = c(rbinom(30, 1, 0.6), rbinom(30, 1, 0.4))
  )

  # Test different alternatives for survival
  result_greater <- TestsSurvBinary(test_data, "OS", "treatment", "control", "greater")
  result_less <- TestsSurvBinary(test_data, "OS", "treatment", "control", "less")
  result_two_sided <- TestsSurvBinary(test_data, "OS", "treatment", "control", "two.sided")

  expect_true(all(c(result_greater$pvalue, result_less$pvalue, result_two_sided$pvalue) >= 0))
  expect_true(all(c(result_greater$pvalue, result_less$pvalue, result_two_sided$pvalue) <= 1))

  # Test different alternatives for binary
  result_or_greater <- TestsSurvBinary(test_data, "OR", "treatment", "control", "greater")
  result_or_less <- TestsSurvBinary(test_data, "OR", "treatment", "control", "less")
  result_or_two_sided <- TestsSurvBinary(test_data, "OR", "treatment", "control", "two.sided")

  expect_true(all(c(result_or_greater$pvalue, result_or_less$pvalue, result_or_two_sided$pvalue) >= 0))
  expect_true(all(c(result_or_greater$pvalue, result_or_less$pvalue, result_or_two_sided$pvalue) <= 1))
})

test_that("TestsSurvBinary handles edge cases", {
  # Empty data
  empty_data <- data.frame(
    ARM = character(0),
    OS = numeric(0),
    Event.OS = numeric(0)
  )

  result_empty <- TestsSurvBinary(empty_data, "OS", "treatment", "control", "greater")
  expect_equal(result_empty$n_treatment, 0)
  expect_equal(result_empty$n_control, 0)
  expect_equal(result_empty$pvalue, 1.0)

  # No events
  no_events_data <- data.frame(
    ARM = c(rep("treatment", 20), rep("control", 20)),
    OS = c(rexp(20, rate = log(2)/15), rexp(20, rate = log(2)/12)),
    Event.OS = 0  # No events
  )

  result_no_events <- TestsSurvBinary(no_events_data, "OS", "treatment", "control", "greater")
  expect_equal(result_no_events$events_treatment, 0)
  expect_equal(result_no_events$events_control, 0)
  expect_equal(result_no_events$pvalue, 1.0)

  # All responses for binary outcome
  all_response_data <- data.frame(
    ARM = c(rep("treatment", 20), rep("control", 20)),
    OR = 1  # All responses
  )

  result_all_response <- TestsSurvBinary(all_response_data, "OR", "treatment", "control", "greater")
  expect_equal(result_all_response$pvalue, 1.0)
})

test_that("TestsSurvBinary handles missing arms", {
  # Data with only one arm
  single_arm_data <- data.frame(
    ARM = rep("treatment", 30),
    OS = rexp(30, rate = log(2)/15),
    Event.OS = 1
  )

  result_single <- TestsSurvBinary(single_arm_data, "OS", "treatment", "control", "greater")
  expect_equal(result_single$n_treatment, 30)
  expect_equal(result_single$n_control, 0)
  expect_equal(result_single$pvalue, 1.0)
})

test_that("TestsSurvBinary calculates event counts correctly", {
  test_data <- data.frame(
    ARM = c(rep("A", 25), rep("B", 25)),
    OS = rexp(50, rate = log(2)/12),
    PFS = rexp(50, rate = log(2)/8),
    OR = c(rep(1, 15), rep(0, 10), rep(1, 10), rep(0, 15)),  # 15 events in A, 10 in B
    Event.OS = c(rep(1, 20), rep(0, 5), rep(1, 18), rep(0, 7)),  # 20 events in A, 18 in B
    Event.PFS = c(rep(1, 22), rep(0, 3), rep(1, 20), rep(0, 5))   # 22 events in A, 20 in B
  )

  # Test event counting for survival outcomes
  result_os <- TestsSurvBinary(test_data, "OS", "A", "B", "greater")
  expect_equal(result_os$events_treatment, 20)
  expect_equal(result_os$events_control, 18)

  result_pfs <- TestsSurvBinary(test_data, "PFS", "A", "B", "greater")
  expect_equal(result_pfs$events_treatment, 22)
  expect_equal(result_pfs$events_control, 20)

  # Test event counting for binary outcome
  result_or <- TestsSurvBinary(test_data, "OR", "A", "B", "greater")
  expect_equal(result_or$events_treatment, 15)
  expect_equal(result_or$events_control, 10)
})
