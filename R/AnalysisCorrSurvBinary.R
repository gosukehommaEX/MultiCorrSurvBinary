#' Advanced Event-Driven Analysis for Multiple Arms and Subgroups
#'
#' This function performs comprehensive event-driven analysis on simulated data from
#' rCorrSurvBinaryMultiArmSubgroup, with subgroup prioritization and multiple arm comparisons.
#'
#' @param data A tibble from rCorrSurvBinaryMultiArmSubgroup containing simulation results
#' @param E A numeric vector of event numbers for sequential analyses (e.g., c(100, 200, 300))
#' @param prioritize Character string specifying which outcome drives the analysis timing.
#'   Must be either "OS" or "PFS" and should match the prioritize parameter used
#'   in data generation.
#' @param subgroup.prioritize Character vector specifying subgroup analysis priority.
#'   Can include 'entire' for overall population and specific subgroup names.
#'   Example: c('entire', 'sub1', 'sub2'). Analysis timing is determined by the
#'   first subgroup in this vector.
#' @param alternative Direction of alternative hypothesis. Options are:
#'   "greater" (default): treatment > control (improvement in treatment)
#'   "less": treatment < control
#'   "two.sided": treatment ≠ control
#'
#' @return A tibble containing analysis results with columns:
#'   \item{sim}{Simulation iteration number}
#'   \item{analysis_event}{Event number triggering the analysis (from E vector)}
#'   \item{timing_subgroup}{Subgroup used for analysis timing (first in subgroup.prioritize)}
#'   \item{analysis_subgroup}{Subgroup being analyzed}
#'   \item{outcome}{Outcome being analyzed (OS, PFS, or OR)}
#'   \item{comparison}{Comparison being made (e.g., "arm1_vs_arm3")}
#'   \item{treatment_arm}{Treatment arm in comparison}
#'   \item{control_arm}{Control arm in comparison}
#'   \item{analysis_time}{Analysis time when E-th event occurred}
#'   \item{n_treatment}{Sample size in treatment arm at analysis time}
#'   \item{n_control}{Sample size in control arm at analysis time}
#'   \item{events_treatment, events_control}{Events/responses in each arm}
#'   \item{pvalue}{P-value from statistical test}
#'
#' @details
#' The algorithm performs the following:
#'
#' 1. **Timing Determination**: Uses first subgroup in subgroup.prioritize to determine analysis timing
#' 2. **Multiple Comparisons**: All treatment arms vs. control arm (highest number)
#' 3. **Subgroup Analysis**: Analyzes each subgroup specified in subgroup.prioritize
#' 4. **Entire Population**: When "entire" is specified, combines all subgroups
#' 5. **Statistical Tests**:
#'    - OS/PFS: One-sided log-rank test
#'    - OR: One-sided chi-squared test with pooled proportion
#'
#' Total number of tests per simulation =
#' length(subgroup.prioritize) × available_outcomes × (n_arms - 1)
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(dplyr)
#' library(tibble)
#' library(survival)
#' library(copula)
#'
#' # Generate sample data
#' arm.params.subgroups <- list(
#'   arm1 = list(
#'     sub1 = list(
#'       mst.OS = 18, mst.PFS = 10, p.OR = 0.5, n = 150,
#'       rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
#'     ),
#'     sub2 = list(
#'       mst.OS = 14, mst.PFS = 7, p.OR = 0.4, n = 50,
#'       rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
#'     )
#'   ),
#'   arm2 = list(
#'     sub1 = list(
#'       mst.OS = 12, mst.PFS = 5, p.OR = 0.3, n = 150,
#'       rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
#'     ),
#'     sub2 = list(
#'       mst.OS = 10, mst.PFS = 5, p.OR = 0.2, n = 50,
#'       rho.OS.PFS = 0.4, rho.OS.OR = 0.2, rho.PFS.OR = 0.3
#'     )
#'   )
#' )
#'
#' result_data <- rCorrSurvBinaryMultiArmSubgroup(
#'   nsim = 1000,
#'   outcomes = c('OS', 'PFS', 'OR'),
#'   arm.params = arm.params.subgroups,
#'   tau = 18,
#'   seed = 456,
#'   prioritize = 'PFS'
#' )
#'
#' # Perform event-driven analysis
#' analysis_results <- AnalysisCorrSurvBinary(
#'   data = result_data,
#'   E = c(100, 200, 300),
#'   prioritize = "PFS",
#'   subgroup.prioritize = c("entire", "sub1"),
#'   alternative = "greater"
#' )
#'
#' # Summary by subgroup and outcome
#' analysis_results %>%
#'   group_by(analysis_event, analysis_subgroup, outcome) %>%
#'   summarise(
#'     n_comparisons = n(),
#'     significant_tests = sum(pvalue < 0.05, na.rm = TRUE),
#'     mean_analysis_time_years = mean(analysis_time, na.rm = TRUE),
#'     mean_analysis_time_months = mean(analysis_time * 12, na.rm = TRUE),
#'     mean_sample_size = mean(n_treatment + n_control, na.rm = TRUE),
#'     .groups = 'drop'
#'   )
#' }
#'
#' @import dplyr
#' @import survival
#' @importFrom stats pchisq pnorm
#' @export
AnalysisCorrSurvBinary <- function(data,
                                   E,
                                   prioritize,
                                   subgroup.prioritize,
                                   alternative = "greater") {

  # Validate inputs
  if (!prioritize %in% c("OS", "PFS")) {
    stop("prioritize must be either 'OS' or 'PFS'")
  }

  if (!alternative %in% c("greater", "less", "two.sided")) {
    stop("alternative must be 'greater', 'less', or 'two.sided'")
  }

  if (length(subgroup.prioritize) == 0) {
    stop("subgroup.prioritize must contain at least one element")
  }

  # Check if prioritized outcome exists in data
  if (!prioritize %in% names(data)) {
    stop(paste("Prioritized outcome", prioritize, "not found in data"))
  }

  # Determine available outcomes and arms
  available_outcomes <- intersect(c("OS", "PFS", "OR"), names(data))
  available_arms <- sort(unique(data$ARM))
  available_subgroups <- sort(unique(data$SUBGROUP[!is.na(data$SUBGROUP)]))

  # Validate subgroup.prioritize
  valid_subgroups <- c("entire", available_subgroups)
  invalid_subgroups <- setdiff(subgroup.prioritize, valid_subgroups)
  if (length(invalid_subgroups) > 0) {
    stop(paste("Invalid subgroups in subgroup.prioritize:", paste(invalid_subgroups, collapse = ", ")))
  }

  # Determine control arm (highest number)
  control_arm <- available_arms[length(available_arms)]
  treatment_arms <- available_arms[-length(available_arms)]

  cat(sprintf("Available arms: %s\n", paste(available_arms, collapse = ", ")))
  cat(sprintf("Control arm: %s\n", control_arm))
  cat(sprintf("Treatment arms: %s\n", paste(treatment_arms, collapse = ", ")))
  cat(sprintf("Available outcomes: %s\n", paste(available_outcomes, collapse = ", ")))
  cat(sprintf("Subgroup prioritization: %s\n", paste(subgroup.prioritize, collapse = ", ")))

  # Calculate total number of tests
  n_comparisons <- length(treatment_arms)
  n_subgroups <- length(subgroup.prioritize)
  n_outcomes <- length(available_outcomes)
  total_tests_per_sim <- n_comparisons * n_subgroups * n_outcomes * length(E)

  cat(sprintf("Total tests per simulation: %d\n", total_tests_per_sim))

  # Initialize results list
  all_results <- list()

  # Get unique simulation numbers
  unique_sims <- unique(data$sim)
  n_sims <- length(unique_sims)

  cat(sprintf("Performing analysis for %d simulations...\n", n_sims))

  # Progress tracking
  progress_interval <- max(1, floor(n_sims / 20))

  # Main analysis loop
  for (sim_idx in seq_along(unique_sims)) {
    current_sim <- unique_sims[sim_idx]

    # Progress reporting
    if (sim_idx %% progress_interval == 0 || sim_idx == n_sims) {
      cat(sprintf("Processing simulation %d/%d (%.1f%%)\n",
                  sim_idx, n_sims, 100 * sim_idx / n_sims))
    }

    # Extract data for current simulation
    sim_data <- data %>%
      filter(sim == current_sim) %>%
      arrange(ARM, SUBGROUP, patientID)

    # Calculate total times
    if ("OS" %in% available_outcomes) {
      sim_data$Total.OS <- sim_data$Accrual + sim_data$OS
    }
    if ("PFS" %in% available_outcomes) {
      sim_data$Total.PFS <- sim_data$Accrual + sim_data$PFS
    }

    # Get timing subgroup (first in subgroup.prioritize)
    timing_subgroup <- subgroup.prioritize[1]

    # Prepare timing data
    if (timing_subgroup == "entire") {
      timing_data <- sim_data
    } else {
      timing_data <- sim_data %>% filter(SUBGROUP == timing_subgroup)
    }

    # Skip if no data for timing subgroup
    if (nrow(timing_data) == 0) {
      warning(sprintf("Simulation %d: No data for timing subgroup %s", current_sim, timing_subgroup))
      next
    }

    # Sort timing data by prioritized total time
    prioritized_total_col <- paste0("Total.", prioritize)
    timing_data <- timing_data %>% arrange(!!sym(prioritized_total_col))

    # Analyze at each event number
    for (e_val in E) {

      # Skip if not enough patients for this event number in timing subgroup
      if (nrow(timing_data) < e_val) {
        warning(sprintf("Simulation %d: Not enough patients in timing subgroup %s (%d) for event number %d",
                        current_sim, timing_subgroup, nrow(timing_data), e_val))
        next
      }

      # Get the analysis time (time of E-th event in timing subgroup)
      analysis_time <- timing_data[[prioritized_total_col]][e_val]

      # For each subgroup to analyze
      for (analysis_subgroup in subgroup.prioritize) {

        # Prepare analysis data for current subgroup
        if (analysis_subgroup == "entire") {
          subgroup_data <- sim_data
        } else {
          subgroup_data <- sim_data %>% filter(SUBGROUP == analysis_subgroup)
        }

        # Skip if no data for analysis subgroup
        if (nrow(subgroup_data) == 0) {
          next
        }

        # Define events and exclude patients not yet enrolled
        subgroup_data <- subgroup_data %>%
          filter(Accrual <= analysis_time) %>%
          mutate(
            Event.OS = if("OS" %in% available_outcomes) ifelse(Total.OS <= analysis_time, 1, 0) else NA,
            Event.PFS = if("PFS" %in% available_outcomes) ifelse(Total.PFS <= analysis_time, 1, 0) else NA
          )

        # For each outcome
        for (outcome in available_outcomes) {

          # For each treatment vs control comparison
          for (treatment_arm in treatment_arms) {

            # Filter data for this comparison
            comparison_data <- subgroup_data %>%
              filter(ARM %in% c(treatment_arm, control_arm))

            # Skip if insufficient data
            if (nrow(comparison_data) == 0) {
              next
            }

            # Perform statistical test
            test_result <- TestsSurvBinary(
              data = comparison_data,
              outcome = outcome,
              treatment_arm = treatment_arm,
              control_arm = control_arm,
              alternative = alternative
            )

            # Create result row
            result_row <- data.frame(
              sim = current_sim,
              analysis_event = e_val,
              timing_subgroup = timing_subgroup,
              analysis_subgroup = analysis_subgroup,
              outcome = outcome,
              comparison = paste(treatment_arm, "vs", control_arm, sep = "_"),
              treatment_arm = treatment_arm,
              control_arm = control_arm,
              analysis_time = analysis_time,
              n_treatment = test_result$n_treatment,
              n_control = test_result$n_control,
              events_treatment = test_result$events_treatment,
              events_control = test_result$events_control,
              pvalue = test_result$pvalue,
              stringsAsFactors = FALSE
            )

            all_results[[length(all_results) + 1]] <- result_row
          }
        }
      }
    }
  }

  # Combine all results
  if (length(all_results) == 0) {
    warning("No results generated")
    return(tibble())
  }

  final_results <- do.call(rbind, all_results)

  cat("Event-driven analysis completed.\n")
  cat(sprintf("Total results generated: %d\n", nrow(final_results)))

  return(dplyr::as_tibble(final_results))
}
