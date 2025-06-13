#' Generate Correlated Time-to-Event and Binary Outcomes for Multiple Arms and Subgroups
#'
#' This function extends rCorrSurvBinary to handle multiple treatment arms and
#' subgroups, generating correlated time-to-event and binary outcomes including
#' overall survival (OS), progression-free survival (PFS), and objective response (OR).
#'
#' @param nsim A numeric value representing number of iterations
#' @param outcomes A character vector including 'OS', 'PFS' and/or 'OR' (same for all arms/subgroups)
#' @param arm.params A named list where each element represents a treatment arm.
#'   Each arm contains either:
#'   - Direct parameters (for single population): mst.OS, mst.PFS, p.OR, rho.OS.PFS, rho.OS.OR, rho.PFS.OR, n
#'   - Named list of subgroups (for multiple subgroups): sub1 = list(...), sub2 = list(...), etc.
#'   Each subgroup list should contain: mst.OS, mst.PFS, p.OR, rho.OS.PFS, rho.OS.OR, rho.PFS.OR, n
#' @param tau A numeric value representing patients' accrual period in years
#'   (patients will be enrolled within [0, tau], default: 12)
#' @param seed A numeric value representing seed number for reproducing results
#' @param validate.bounds Logical, whether to validate correlation bounds against
#'   Fréchet-Hoeffding bounds before simulation (default: TRUE)
#' @param prioritize Character string specifying which distribution to preserve
#'   when applying OS >= PFS constraint. Options are "OS" (default) or "PFS".
#'
#' @return A tibble containing simulation results with columns:
#'   \item{sim}{Simulation iteration number (1 to nsim)}
#'   \item{ARM}{Treatment arm identifier (arm1, arm2, ...)}
#'   \item{SUBGROUP}{Subgroup identifier (sub1, sub2, ... or NA if no subgroups)}
#'   \item{patientID}{Patient identifier within each arm/subgroup/iteration}
#'   \item{OS}{Overall survival time in years (if 'OS' in outcomes)}
#'   \item{PFS}{Progression-free survival time in years (if 'PFS' in outcomes)}
#'   \item{OR}{Objective response indicator 0/1 (if 'OR' in outcomes)}
#'   \item{Accrual}{Patient accrual time in years}
#'
#' @details
#' The function supports two main scenarios:
#' \itemize{
#'   \item Single population per arm: Each arm has one set of parameters
#'   \item Multiple subgroups per arm: Each arm contains multiple subgroups with different parameters
#' }
#'
#' Parameter structure examples:
#' \itemize{
#'   \item Single population: arm.params = list(arm1 = list(mst.OS = 12, mst.PFS = 6, p.OR = 0.4, n = 100, ...))
#'   \item Multiple subgroups: arm.params = list(arm1 = list(sub1 = list(mst.OS = 12, n = 50, ...), sub2 = list(mst.OS = 15, n = 50, ...)))
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Two arms without subgroups
#' arm.params.simple <- list(
#'   arm1 = list(
#'     mst.OS = 18, mst.PFS = 12, p.OR = 0.6, n = 100,
#'     rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
#'   ),
#'   arm2 = list(
#'     mst.OS = 12, mst.PFS = 8, p.OR = 0.4, n = 100,
#'     rho.OS.PFS = 0.5, rho.OS.OR = 0.3, rho.PFS.OR = 0.4
#'   )
#' )
#'
#' result1 <- rCorrSurvBinaryMultiArmSubgroup(
#'   nsim = 1000,
#'   outcomes = c('OS', 'PFS', 'OR'),
#'   arm.params = arm.params.simple,
#'   tau = 24,
#'   seed = 123
#' )
#'
#' # Example 2: Two arms with two subgroups each
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
#' result2 <- rCorrSurvBinaryMultiArmSubgroup(
#'   nsim = 500,
#'   outcomes = c('OS', 'PFS', 'OR'),
#'   arm.params = arm.params.subgroups,
#'   tau = 18,
#'   seed = 456
#' )
#' }
#'
#' @import copula
#' @import dplyr
#' @importFrom stats qexp qbinom runif rexp rbinom
#' @export
rCorrSurvBinaryMultiArmSubgroup <- function(nsim,
                                            outcomes,
                                            arm.params,
                                            tau = 12,
                                            seed = NULL,
                                            validate.bounds = TRUE,
                                            prioritize = "OS") {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (!is.list(arm.params) || is.null(names(arm.params))) {
    stop("arm.params must be a named list")
  }

  # Validate outcomes
  valid.outcomes <- c("OS", "PFS", "OR")
  if (!all(outcomes %in% valid.outcomes)) {
    stop("outcomes must contain only 'OS', 'PFS', and/or 'OR'")
  }

  # Validate prioritize parameter
  if (!prioritize %in% c("OS", "PFS")) {
    stop("prioritize must be either 'OS' or 'PFS'")
  }

  # Function to check if a list contains subgroups or direct parameters
  has.subgroups <- function(arm.param) {
    # Check if all elements are lists and have names starting with "sub"
    if (!all(sapply(arm.param, is.list))) return(FALSE)
    subgroup.names <- names(arm.param)
    return(all(grepl("^sub", subgroup.names)))
  }

  # Function to validate parameter completeness
  validate.params <- function(params, context) {
    required.params <- c("n")

    if ("OS" %in% outcomes) required.params <- c(required.params, "mst.OS")
    if ("PFS" %in% outcomes) required.params <- c(required.params, "mst.PFS")
    if ("OR" %in% outcomes) required.params <- c(required.params, "p.OR")

    # Correlation parameters (only check if multiple outcomes)
    if (length(outcomes) > 1) {
      if (all(c("OS", "PFS") %in% outcomes)) required.params <- c(required.params, "rho.OS.PFS")
      if (all(c("OS", "OR") %in% outcomes)) required.params <- c(required.params, "rho.OS.OR")
      if (all(c("PFS", "OR") %in% outcomes)) required.params <- c(required.params, "rho.PFS.OR")
    }

    missing.params <- setdiff(required.params, names(params))
    if (length(missing.params) > 0) {
      stop(sprintf("Missing required parameters in %s: %s",
                   context, paste(missing.params, collapse = ", ")))
    }
  }

  # Parse arm.params structure and validate
  simulation.configs <- list()

  for (arm.name in names(arm.params)) {
    arm.param <- arm.params[[arm.name]]

    if (has.subgroups(arm.param)) {
      # Multiple subgroups case
      for (subgroup.name in names(arm.param)) {
        subgroup.param <- arm.param[[subgroup.name]]
        context <- sprintf("%s -> %s", arm.name, subgroup.name)
        validate.params(subgroup.param, context)

        config <- list(
          arm = arm.name,
          subgroup = subgroup.name,
          params = subgroup.param
        )
        simulation.configs[[length(simulation.configs) + 1]] <- config
      }
    } else {
      # Single population case
      context <- arm.name
      validate.params(arm.param, context)

      config <- list(
        arm = arm.name,
        subgroup = NA,
        params = arm.param
      )
      simulation.configs[[length(simulation.configs) + 1]] <- config
    }
  }

  # Validate correlation bounds for all configurations if requested
  if (validate.bounds) {
    cat("Validating correlation bounds for all arm/subgroup combinations...\n")

    for (i in seq_along(simulation.configs)) {
      config <- simulation.configs[[i]]
      params <- config$params
      context <- if (is.na(config$subgroup)) config$arm else paste(config$arm, config$subgroup, sep = " -> ")

      # Convert parameters for validation
      lambda.OS <- if ("OS" %in% outcomes) log(2) / params$mst.OS else NULL
      lambda.PFS <- if ("PFS" %in% outcomes) log(2) / params$mst.PFS else NULL
      p.OR <- if ("OR" %in% outcomes) params$p.OR else NULL

      rho.OS.PFS <- if ("rho.OS.PFS" %in% names(params)) params$rho.OS.PFS else 0
      rho.OS.OR <- if ("rho.OS.OR" %in% names(params)) params$rho.OS.OR else 0
      rho.PFS.OR <- if ("rho.PFS.OR" %in% names(params)) params$rho.PFS.OR else 0

      validation <- CorrBounds(
        outcomes = outcomes,
        lambda.OS = lambda.OS,
        lambda.PFS = lambda.PFS,
        p.OR = p.OR,
        rho.OS.PFS = rho.OS.PFS,
        rho.OS.OR = rho.OS.OR,
        rho.PFS.OR = rho.PFS.OR
      )

      if (!validation$valid) {
        cat(sprintf("\n=== CORRELATION VALIDATION FAILED FOR %s ===\n", context))
        for (error in validation$errors) {
          cat("ERROR:", error, "\n")
        }

        cat("\n=== TRUE FRÉCHET-HOEFFDING BOUNDS ===\n")
        for (bound.name in names(validation$bounds)) {
          bound.info <- validation$bounds[[bound.name]]
          cat(sprintf("%s: [%.6f, %.6f]\n",
                      bound.name, bound.info$lower, bound.info$upper))
        }

        stop(sprintf("Correlation validation failed for %s. Please adjust correlation parameters.", context))
      }
    }

    cat("All correlation bounds validation passed.\n")
  }

  # Main simulation loop
  result.list <- list()

  for (sim.iter in 1:nsim) {
    sim.results <- list()

    for (config in simulation.configs) {
      params <- config$params

      # Generate data for this arm/subgroup combination
      sim.data <- rCorrSurvBinary(
        nsim = 1,
        outcomes = outcomes,
        n = params$n,
        mst.OS = if ("OS" %in% outcomes) params$mst.OS else NULL,
        mst.PFS = if ("PFS" %in% outcomes) params$mst.PFS else NULL,
        p.OR = if ("OR" %in% outcomes) params$p.OR else NULL,
        rho.OS.PFS = if ("rho.OS.PFS" %in% names(params)) params$rho.OS.PFS else 0,
        rho.OS.OR = if ("rho.OS.OR" %in% names(params)) params$rho.OS.OR else 0,
        rho.PFS.OR = if ("rho.PFS.OR" %in% names(params)) params$rho.PFS.OR else 0,
        tau = tau,
        seed = NULL,  # Don't set seed for individual calls to maintain randomness
        validate.bounds = FALSE,  # Already validated above
        prioritize = prioritize
      )

      # Add ARM and SUBGROUP columns
      sim.data$ARM <- config$arm
      sim.data$SUBGROUP <- config$subgroup

      # Update sim number
      sim.data$sim <- sim.iter

      # Reorder columns
      col.order <- c("sim", "ARM", "SUBGROUP", "patientID",
                     intersect(c("OS", "PFS", "OR"), names(sim.data)),
                     "Accrual")
      sim.data <- sim.data[, col.order]

      sim.results[[length(sim.results) + 1]] <- sim.data
    }

    # Combine results for this simulation
    result.list[[sim.iter]] <- do.call(rbind, sim.results)
  }

  # Combine all simulation results
  final.result <- do.call(rbind, result.list)

  # Reassign patient IDs to be unique within each simulation, arm, and subgroup
  final.result <- final.result %>%
    group_by(sim, ARM, SUBGROUP) %>%
    mutate(patientID = row_number()) %>%
    ungroup()

  # Convert to tibble and return
  return(dplyr::as_tibble(final.result))
}
