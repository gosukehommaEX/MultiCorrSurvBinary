#' Generate Correlated Time-to-Event and Binary Outcomes
#'
#' This function generates pseudo random numbers of correlated time-to-event
#' and binary outcomes including overall survival (OS), progression-free
#' survival (PFS), and objective response (OR). The function uses a normal
#' copula to model correlation structures while maintaining specified marginal
#' distributions.
#'
#' @param nsim A numeric value representing number of iterations
#' @param outcomes A character vector including 'OS', 'PFS' and/or 'OR'
#' @param n A numeric value representing sample size per iteration
#' @param mst.OS A numeric value representing median survival time of OS in years
#' @param mst.PFS A numeric value representing median survival time of PFS in years
#' @param p.OR A numeric value representing probability of OR (0 <= p.OR <= 1)
#' @param rho.OS.PFS A numeric value representing correlation between OS and PFS
#'   (must be within Frechet-Hoeffding bounds, default: 0)
#' @param rho.OS.OR A numeric value representing correlation between OS and OR
#'   (must be within Frechet-Hoeffding bounds, default: 0)
#' @param rho.PFS.OR A numeric value representing correlation between PFS and OR
#'   (must be within Frechet-Hoeffding bounds, default: 0)
#' @param tau A numeric value representing patients' accrual period in years
#'   (patients will be enrolled within [0, tau], default: 12)
#' @param seed A numeric value representing seed number for reproducing results
#' @param validate.bounds Logical, whether to validate correlation bounds against
#'   Frechet-Hoeffding bounds before simulation (default: TRUE)
#' @param prioritize Character string specifying which distribution to preserve
#'   when applying OS >= PFS constraint. Options are "OS" (default) or "PFS".
#'   Only applies when both OS and PFS are included in outcomes.
#'
#' @return A tibble containing simulation results with columns:
#'   \item{sim}{Simulation iteration number (1 to nsim)}
#'   \item{patientID}{Patient identifier within each iteration (1 to n)}
#'   \item{OS}{Overall survival time in years (if 'OS' in outcomes)}
#'   \item{PFS}{Progression-free survival time in years (if 'PFS' in outcomes)}
#'   \item{OR}{Objective response indicator 0/1 (if 'OR' in outcomes)}
#'   \item{Accrual}{Patient accrual time in years}
#'
#' @details
#' The function implements the following assumptions:
#' \itemize{
#'   \item Marginal distribution of OS follows exponential distribution
#'   \item Marginal distribution of PFS follows exponential distribution
#'   \item Marginal distribution of OR follows Bernoulli distribution
#'   \item OS >= PFS constraint is enforced when both outcomes are included
#'   \item Patient accrual follows uniform distribution over [0, tau]
#'   \item Correlation structure is modeled using normal copula
#' }
#'
#' When both OS and PFS are included, the constraint OS >= PFS is enforced.
#' The prioritize parameter determines which distribution is preserved:
#' \itemize{
#'   \item prioritize = "OS": PFS is adjusted to ensure PFS <= OS (OS distribution preserved)
#'   \item prioritize = "PFS": OS is adjusted to ensure OS >= PFS (PFS distribution preserved)
#' }
#'
#' @examples
#' # Generate correlated OS, PFS, and OR
#' result <- rCorrSurvBinary(
#'   nsim = 100,
#'   outcomes = c('OS', 'PFS', 'OR'),
#'   n = 200,
#'   mst.OS = 12,
#'   mst.PFS = 6,
#'   p.OR = 0.4,
#'   rho.OS.PFS = 0.5,
#'   rho.OS.OR = 0.3,
#'   rho.PFS.OR = 0.4,
#'   tau = 24,
#'   seed = 123,
#'   prioritize = "OS"
#' )
#'
#' # Generate only OS and PFS with correlation, preserving PFS distribution
#' result2 <- rCorrSurvBinary(
#'   nsim = 100,
#'   outcomes = c('OS', 'PFS'),
#'   n = 100,
#'   mst.OS = 15,
#'   mst.PFS = 8,
#'   rho.OS.PFS = 0.7,
#'   tau = 18,
#'   seed = 456,
#'   prioritize = "PFS"
#' )
#'
#' @import copula
#' @import dplyr
#' @importFrom stats qexp qbinom runif rexp rbinom
#' @export
rCorrSurvBinary <- function(nsim,
                            outcomes,
                            n,
                            mst.OS = NULL,
                            mst.PFS = NULL,
                            p.OR = NULL,
                            rho.OS.PFS = 0,
                            rho.OS.OR = 0,
                            rho.PFS.OR = 0,
                            tau = 12,
                            seed = NULL,
                            validate.bounds = TRUE,
                            prioritize = "OS") {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate prioritize parameter
  if (!prioritize %in% c("OS", "PFS")) {
    stop("prioritize must be either 'OS' or 'PFS'")
  }

  # Check if prioritize is relevant (only when both OS and PFS are present)
  if (!all(c("OS", "PFS") %in% outcomes) && prioritize != "OS") {
    warning("prioritize parameter is only relevant when both OS and PFS are in outcomes. Using default behavior.")
  }

  # Validate outcomes
  valid.outcomes <- c("OS", "PFS", "OR")
  if (!all(outcomes %in% valid.outcomes)) {
    stop("outcomes must contain only 'OS', 'PFS', and/or 'OR'")
  }

  # Validate required parameters based on outcomes
  if ("OS" %in% outcomes && is.null(mst.OS)) {
    stop("mst.OS must be specified when 'OS' is in outcomes")
  }
  if ("PFS" %in% outcomes && is.null(mst.PFS)) {
    stop("mst.PFS must be specified when 'PFS' is in outcomes")
  }
  if ("OR" %in% outcomes && is.null(p.OR)) {
    stop("p.OR must be specified when 'OR' is in outcomes")
  }

  # Convert median survival times to hazard rates
  lambda.OS <- if ("OS" %in% outcomes) log(2) / mst.OS else NULL
  lambda.PFS <- if ("PFS" %in% outcomes) log(2) / mst.PFS else NULL

  # Validate correlation bounds if requested
  if (validate.bounds) {
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
      # Print detailed error information with bounds
      cat("=== CORRELATION VALIDATION FAILED ===\n")
      for (error in validation$errors) {
        cat("ERROR:", error, "\n")
      }

      cat("\n=== TRUE FRECHET-HOEFFDING BOUNDS ===\n")
      for (bound.name in names(validation$bounds)) {
        bound.info <- validation$bounds[[bound.name]]
        cat(sprintf("%s: [%.6f, %.6f] (%s)\n",
                    bound.name, bound.info$lower, bound.info$upper, bound.info$method))
      }

      cat(sprintf("\nBounds computation time: %.3f seconds\n", validation$computation.time))

      stop("Correlation validation failed. Please adjust correlation parameters to be within true Frechet-Hoeffding bounds.")
    }
  }

  # Initialize result list
  result.list <- list()

  # Main simulation loop
  for (sim.iter in 1:nsim) {

    # Generate accrual times
    accrual.times <- runif(n, min = 0, max = tau)

    # Generate correlated outcomes based on the number of outcomes
    n.outcomes <- length(outcomes)

    if (n.outcomes == 1) {
      # Single outcome case
      if (outcomes == "OS") {
        os.values <- rexp(n, rate = lambda.OS)
        sim.data <- data.frame(
          sim = sim.iter,
          patientID = 1:n,
          OS = os.values,
          Accrual = accrual.times
        )
      } else if (outcomes == "PFS") {
        pfs.values <- rexp(n, rate = lambda.PFS)
        sim.data <- data.frame(
          sim = sim.iter,
          patientID = 1:n,
          PFS = pfs.values,
          Accrual = accrual.times
        )
      } else if (outcomes == "OR") {
        or.values <- rbinom(n, size = 1, prob = p.OR)
        sim.data <- data.frame(
          sim = sim.iter,
          patientID = 1:n,
          OR = or.values,
          Accrual = accrual.times
        )
      }

    } else {
      # Multiple outcomes case - use copula for correlation

      # Create correlation matrix
      if (n.outcomes == 2) {
        if (all(c("OS", "PFS") %in% outcomes)) {
          corr.matrix <- matrix(c(1, rho.OS.PFS, rho.OS.PFS, 1), nrow = 2)
        } else if (all(c("OS", "OR") %in% outcomes)) {
          corr.matrix <- matrix(c(1, rho.OS.OR, rho.OS.OR, 1), nrow = 2)
        } else if (all(c("PFS", "OR") %in% outcomes)) {
          corr.matrix <- matrix(c(1, rho.PFS.OR, rho.PFS.OR, 1), nrow = 2)
        }
      } else if (n.outcomes == 3) {
        corr.matrix <- matrix(c(1, rho.OS.PFS, rho.OS.OR,
                                rho.OS.PFS, 1, rho.PFS.OR,
                                rho.OS.OR, rho.PFS.OR, 1),
                              nrow = 3, ncol = 3)
      }

      # Check if correlation matrix is positive definite
      if (!all(eigen(corr.matrix)$values > 0)) {
        stop("The specified correlation coefficients do not form a positive definite matrix")
      }

      # Generate correlated uniform random variables using normal copula
      normal.cop <- normalCopula(dim = n.outcomes, dispstr = "un")
      normal.cop@parameters <- corr.matrix[upper.tri(corr.matrix)]
      uniform.random <- rCopula(n, normal.cop)

      # Initialize data frame
      sim.data <- data.frame(
        sim = sim.iter,
        patientID = 1:n,
        Accrual = accrual.times
      )

      # Transform to target distributions
      col.index <- 1

      if ("OS" %in% outcomes) {
        os.values <- qexp(uniform.random[, col.index], rate = lambda.OS)
        sim.data$OS <- os.values
        col.index <- col.index + 1
      }

      if ("PFS" %in% outcomes) {
        pfs.initial <- qexp(uniform.random[, col.index], rate = lambda.PFS)
        sim.data$PFS <- pfs.initial
        col.index <- col.index + 1
      }

      if ("OR" %in% outcomes) {
        or.values <- qbinom(uniform.random[, col.index], size = 1, prob = p.OR)
        sim.data$OR <- or.values
      }

      # Apply OS >= PFS constraint if both are present
      if (all(c("OS", "PFS") %in% outcomes)) {
        if (prioritize == "OS") {
          # Preserve OS distribution, adjust PFS
          sim.data$PFS <- pmin(sim.data$PFS, sim.data$OS)
        } else if (prioritize == "PFS") {
          # Preserve PFS distribution, adjust OS
          sim.data$OS <- pmax(sim.data$OS, sim.data$PFS)
        }
      }
    }

    result.list[[sim.iter]] <- sim.data
  }

  # Combine all simulation results
  final.result <- do.call(rbind, result.list)

  # Convert to tibble and return
  return(dplyr::as_tibble(final.result))
}
