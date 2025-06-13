#' Perform Statistical Test for Single Outcome and Comparison
#'
#' This function performs appropriate statistical tests for different outcome types
#' including time-to-event outcomes (OS, PFS) using log-rank tests and binary
#' outcomes (OR) using chi-squared tests with pooled proportion.
#'
#' @param data Comparison dataset (two arms only)
#' @param outcome Outcome name ("OS", "PFS", or "OR")
#' @param treatment_arm Treatment arm name
#' @param control_arm Control arm name
#' @param alternative Alternative hypothesis direction ("greater", "less", "two.sided")
#' @return List with test results including sample sizes, events, and p-value
#'
#' @details
#' For time-to-event outcomes (OS, PFS):
#' - Uses survdiff() from survival package for log-rank test
#' - Calculates one-sided p-values by comparing observed vs expected events
#' - "greater" alternative tests if treatment arm has better survival (fewer events)
#'
#' For binary outcomes (OR):
#' - Uses pooled proportion chi-squared test
#' - Calculates z-statistic with pooled standard error
#' - "greater" alternative tests if treatment arm has higher response rate
#'
#' @import survival
#' @importFrom stats pchisq pnorm
#' @export
TestsSurvBinary <- function(data, outcome, treatment_arm, control_arm, alternative) {

  # Get sample sizes
  n_treatment <- sum(data$ARM == treatment_arm)
  n_control <- sum(data$ARM == control_arm)

  # Initialize result
  result <- list(
    n_treatment = n_treatment,
    n_control = n_control,
    events_treatment = 0,
    events_control = 0,
    pvalue = 1.0
  )

  # Skip if no patients in either arm
  if (n_treatment == 0 || n_control == 0) {
    return(result)
  }

  if (outcome %in% c("OS", "PFS")) {
    # Time-to-event analysis using log-rank test
    event_col <- paste0("Event.", outcome)

    # Get event counts
    result$events_treatment <- sum(data[data$ARM == treatment_arm, event_col], na.rm = TRUE)
    result$events_control <- sum(data[data$ARM == control_arm, event_col], na.rm = TRUE)

    # Skip if no events
    if (result$events_treatment + result$events_control == 0) {
      return(result)
    }

    # Perform log-rank test
    tryCatch({
      surv_obj <- Surv(time = data[[outcome]], event = data[[event_col]])
      surv_fit <- survdiff(surv_obj ~ ARM, data = data)

      chi_stat <- surv_fit$chisq

      # Get observed and expected for treatment arm
      # survdiff orders by ARM alphabetically, so we need to find the right index
      arm_levels <- sort(unique(data$ARM))
      treatment_idx <- which(arm_levels == treatment_arm)

      observed_treatment <- surv_fit$obs[treatment_idx]
      expected_treatment <- surv_fit$exp[treatment_idx]

      # Calculate one-sided p-value based on direction
      if (alternative == "greater") {
        # Treatment better than control (fewer events than expected)
        if (observed_treatment < expected_treatment) {
          result$pvalue <- pchisq(chi_stat, df = 1, lower.tail = FALSE) / 2
        } else {
          result$pvalue <- 1 - pchisq(chi_stat, df = 1, lower.tail = FALSE) / 2
        }
      } else if (alternative == "less") {
        # Treatment worse than control (more events than expected)
        if (observed_treatment > expected_treatment) {
          result$pvalue <- pchisq(chi_stat, df = 1, lower.tail = FALSE) / 2
        } else {
          result$pvalue <- 1 - pchisq(chi_stat, df = 1, lower.tail = FALSE) / 2
        }
      } else {
        # Two-sided test
        result$pvalue <- pchisq(chi_stat, df = 1, lower.tail = FALSE)
      }

    }, error = function(e) {
      warning(paste("Log-rank test failed for", outcome, ":", e$message))
    })

  } else if (outcome == "OR") {
    # Binary outcome analysis using chi-squared test with pooled proportion
    treatment_data <- data[data$ARM == treatment_arm, ]
    control_data <- data[data$ARM == control_arm, ]

    result$events_treatment <- sum(treatment_data$OR, na.rm = TRUE)
    result$events_control <- sum(control_data$OR, na.rm = TRUE)

    total_treatment <- nrow(treatment_data)
    total_control <- nrow(control_data)

    # Perform chi-squared test with pooled proportion
    tryCatch({
      # Pooled proportion under null hypothesis
      pooled_prop <- (result$events_treatment + result$events_control) / (total_treatment + total_control)

      # Avoid edge cases
      if (pooled_prop == 0 || pooled_prop == 1 || total_treatment == 0 || total_control == 0) {
        result$pvalue <- 1.0
      } else {
        # Sample proportions
        prop_treatment <- result$events_treatment / total_treatment
        prop_control <- result$events_control / total_control

        # Standard error using pooled proportion
        se_pooled <- sqrt(pooled_prop * (1 - pooled_prop) * (1/total_treatment + 1/total_control))

        # Z-statistic for difference in proportions
        z_stat <- (prop_treatment - prop_control) / se_pooled

        # Calculate p-value based on alternative hypothesis
        if (alternative == "greater") {
          # Treatment has higher response rate than control
          result$pvalue <- pnorm(z_stat, lower.tail = FALSE)
        } else if (alternative == "less") {
          # Treatment has lower response rate than control
          result$pvalue <- pnorm(z_stat, lower.tail = TRUE)
        } else {
          # Two-sided test
          result$pvalue <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
        }
      }

    }, error = function(e) {
      warning(paste("Chi-squared test failed for OR:", e$message))
    })
  }

  return(result)
}
