#' Calculate and Validate True Fréchet-Hoeffding Correlation Bounds
#'
#' This function calculates the true Fréchet-Hoeffding bounds using analytical
#' integration for different types of variable combinations and validates
#' whether specified correlation coefficients are within these bounds.
#'
#' @param outcomes Character vector of outcomes ('OS', 'PFS', 'OR')
#' @param lambda.OS Rate parameter for OS exponential distribution
#'   (calculated as log(2)/median)
#' @param lambda.PFS Rate parameter for PFS exponential distribution
#'   (calculated as log(2)/median)
#' @param p.OR Probability parameter for OR Bernoulli distribution (0 < p.OR < 1)
#' @param rho.OS.PFS Correlation coefficient between OS and PFS (default: 0)
#' @param rho.OS.OR Correlation coefficient between OS and OR (default: 0)
#' @param rho.PFS.OR Correlation coefficient between PFS and OR (default: 0)
#'
#' @return A list containing:
#'   \item{valid}{Logical indicating if all correlations are valid}
#'   \item{errors}{Character vector of error messages (empty if valid)}
#'   \item{bounds}{List of true Fréchet-Hoeffding bounds for each correlation pair}
#'   \item{eigenvalues}{Eigenvalues of correlation matrix (if applicable)}
#'   \item{computation.time}{Time taken for bounds calculation}
#'
#' @details
#' The function calculates true Fréchet-Hoeffding bounds using analytical integration:
#'
#' \strong{Continuous-Continuous (OS-PFS):}
#' For exponential distributions, the bounds are calculated by finding the
#' correlation of the Fréchet-Hoeffding copulas W (lower bound) and M (upper bound)
#' applied to exponential marginals using R's integrate() function.
#'
#' \strong{Continuous-Binary (OS-OR, PFS-OR):}
#' The bounds are calculated using analytical integration:
#' \deqn{\rho_{max} = \frac{E[XY_{max}] - E[X]E[Y]}{\sigma_X \sigma_Y}}
#' \deqn{\rho_{min} = \frac{E[XY_{min}] - E[X]E[Y]}{\sigma_X \sigma_Y}}
#' where \eqn{Y_{max}} corresponds to binary=1 for largest continuous values,
#' and \eqn{Y_{min}} corresponds to binary=1 for smallest continuous values.
#'
#' \strong{Integration Method:}
#' All bounds are computed using R's integrate() function for high precision
#' analytical integration. No grid approximation is used.
#'
#' @examples
#' # Calculate true bounds for OS, PFS, and OR
#' validation <- CorrBounds(
#'   outcomes = c('OS', 'PFS', 'OR'),
#'   lambda.OS = log(2)/12,  # 12-year median OS
#'   lambda.PFS = log(2)/6,  # 6-year median PFS
#'   p.OR = 0.4,             # 40% response rate
#'   rho.OS.PFS = 0.6,
#'   rho.OS.OR = 0.3,
#'   rho.PFS.OR = 0.4
#' )
#'
#' if (!validation$valid) {
#'   cat("Errors found:\n")
#'   for (error in validation$errors) {
#'     cat("-", error, "\n")
#'   }
#'   cat("\nTrue Fréchet-Hoeffding bounds:\n")
#'   for (bound.name in names(validation$bounds)) {
#'     bound <- validation$bounds[[bound.name]]
#'     cat(sprintf("%s: [%.6f, %.6f]\n", bound.name, bound$lower, bound$upper))
#'   }
#' }
#'
#' @importFrom stats integrate qexp
#' @export
CorrBounds <- function(outcomes,
                       lambda.OS = NULL,
                       lambda.PFS = NULL,
                       p.OR = NULL,
                       rho.OS.PFS = 0,
                       rho.OS.OR = 0,
                       rho.PFS.OR = 0) {

  start.time <- Sys.time()

  # Initialize validation results
  validation.results <- list(
    valid = TRUE,
    errors = character(0),
    bounds = list(),
    eigenvalues = NULL,
    computation.time = NULL
  )

  # Helper function to calculate true bounds for continuous-continuous variables
  calculate.true.continuous.continuous.bounds <- function(lambda1, lambda2) {

    # For exponential distributions, we can calculate Fréchet-Hoeffding bounds
    # using the analytical properties

    # Mean and standard deviation for both distributions
    mean1 <- 1 / lambda1
    mean2 <- 1 / lambda2
    sd1 <- 1 / lambda1
    sd2 <- 1 / lambda2

    # Upper bound: Perfect positive dependence
    # For exponential distributions, this is achieved when both variables
    # use the same uniform random variable u: X1 = F1^(-1)(u), X2 = F2^(-1)(u)

    # E[X1 * X2] under perfect positive dependence
    # = E[qexp(U, λ1) * qexp(U, λ2)] where U ~ Uniform(0,1)

    integrand.upper <- function(u) {
      x1 <- qexp(u, rate = lambda1)
      x2 <- qexp(u, rate = lambda2)
      return(x1 * x2)
    }

    e.x1x2.upper <- integrate(integrand.upper, lower = 0.001, upper = 0.999)$value
    cov.upper <- e.x1x2.upper - mean1 * mean2
    cor.upper <- cov.upper / (sd1 * sd2)

    # Lower bound: Perfect negative dependence
    # X1 = F1^(-1)(U), X2 = F2^(-1)(1-U)

    integrand.lower <- function(u) {
      x1 <- qexp(u, rate = lambda1)
      x2 <- qexp(1 - u + 1e-10, rate = lambda2)  # Small offset to avoid boundary issues
      return(x1 * x2)
    }

    e.x1x2.lower <- integrate(integrand.lower, lower = 0.001, upper = 0.999)$value
    cov.lower <- e.x1x2.lower - mean1 * mean2
    cor.lower <- cov.lower / (sd1 * sd2)

    list(
      type = "Continuous-Continuous (Exponential)",
      lower = cor.lower,
      upper = cor.upper,
      lambda1 = lambda1,
      lambda2 = lambda2,
      method = "Analytical integration using R integrate() function"
    )
  }

  # Helper function to calculate true bounds for continuous-binary variables
  calculate.true.continuous.binary.bounds <- function(lambda.continuous, p.binary) {

    if (p.binary <= 0 || p.binary >= 1) {
      stop("p.binary must be strictly between 0 and 1")
    }

    # Calculate moments of exponential distribution
    mean.continuous <- 1 / lambda.continuous
    sd.continuous <- 1 / lambda.continuous

    # Binary variable moments
    mean.binary <- p.binary
    sd.binary <- sqrt(p.binary * (1 - p.binary))

    # Upper bound: Binary = 1 for largest continuous values
    # Threshold for upper bound: (1-p.binary) quantile
    threshold.upper <- qexp(1 - p.binary, rate = lambda.continuous)

    # E[X * Y] for upper bound: Y=1 when X >= threshold.upper
    integrand.upper <- function(x) {
      x * lambda.continuous * exp(-lambda.continuous * x)
    }

    e.xy.upper <- integrate(integrand.upper,
                            lower = threshold.upper,
                            upper = Inf)$value

    # Covariance and correlation for upper bound
    cov.upper <- e.xy.upper - mean.continuous * mean.binary
    cor.upper <- cov.upper / (sd.continuous * sd.binary)

    # Lower bound: Binary = 1 for smallest continuous values
    threshold.lower <- qexp(p.binary, rate = lambda.continuous)

    # E[X * Y] for lower bound: Y=1 when X <= threshold.lower
    e.xy.lower <- integrate(integrand.upper,
                            lower = 0,
                            upper = threshold.lower)$value

    # Covariance and correlation for lower bound
    cov.lower <- e.xy.lower - mean.continuous * mean.binary
    cor.lower <- cov.lower / (sd.continuous * sd.binary)

    list(
      type = "Continuous-Binary (Exponential-Bernoulli)",
      lower = cor.lower,
      upper = cor.upper,
      p.binary = p.binary,
      lambda.continuous = lambda.continuous,
      threshold.upper = threshold.upper,
      threshold.lower = threshold.lower,
      method = "Analytical integration using R integrate() function"
    )
  }

  # Validate OS-PFS correlation bounds
  if (all(c("OS", "PFS") %in% outcomes)) {
    if (is.null(lambda.OS) || is.null(lambda.PFS)) {
      validation.results$errors <- c(validation.results$errors,
                                     "lambda.OS and lambda.PFS must be specified for OS-PFS correlation")
      validation.results$valid <- FALSE
    } else {
      bounds.OS.PFS <- calculate.true.continuous.continuous.bounds(lambda.OS, lambda.PFS)
      validation.results$bounds$OS.PFS <- bounds.OS.PFS

      if (rho.OS.PFS < bounds.OS.PFS$lower || rho.OS.PFS > bounds.OS.PFS$upper) {
        validation.results$errors <- c(validation.results$errors,
                                       sprintf("rho.OS.PFS (%.6f) is outside true Fréchet-Hoeffding bounds [%.6f, %.6f]",
                                               rho.OS.PFS, bounds.OS.PFS$lower, bounds.OS.PFS$upper))
        validation.results$valid <- FALSE
      }
    }
  }

  # Validate OS-OR correlation bounds
  if (all(c("OS", "OR") %in% outcomes)) {
    if (is.null(lambda.OS) || is.null(p.OR)) {
      validation.results$errors <- c(validation.results$errors,
                                     "lambda.OS and p.OR must be specified for OS-OR correlation")
      validation.results$valid <- FALSE
    } else {
      bounds.OS.OR <- calculate.true.continuous.binary.bounds(lambda.OS, p.OR)
      validation.results$bounds$OS.OR <- bounds.OS.OR

      if (rho.OS.OR < bounds.OS.OR$lower || rho.OS.OR > bounds.OS.OR$upper) {
        validation.results$errors <- c(validation.results$errors,
                                       sprintf("rho.OS.OR (%.6f) is outside true Fréchet-Hoeffding bounds [%.6f, %.6f]",
                                               rho.OS.OR, bounds.OS.OR$lower, bounds.OS.OR$upper))
        validation.results$valid <- FALSE
      }
    }
  }

  # Validate PFS-OR correlation bounds
  if (all(c("PFS", "OR") %in% outcomes)) {
    if (is.null(lambda.PFS) || is.null(p.OR)) {
      validation.results$errors <- c(validation.results$errors,
                                     "lambda.PFS and p.OR must be specified for PFS-OR correlation")
      validation.results$valid <- FALSE
    } else {
      bounds.PFS.OR <- calculate.true.continuous.binary.bounds(lambda.PFS, p.OR)
      validation.results$bounds$PFS.OR <- bounds.PFS.OR

      if (rho.PFS.OR < bounds.PFS.OR$lower || rho.PFS.OR > bounds.PFS.OR$upper) {
        validation.results$errors <- c(validation.results$errors,
                                       sprintf("rho.PFS.OR (%.6f) is outside true Fréchet-Hoeffding bounds [%.6f, %.6f]",
                                               rho.PFS.OR, bounds.PFS.OR$lower, bounds.PFS.OR$upper))
        validation.results$valid <- FALSE
      }
    }
  }

  # Check positive definiteness of correlation matrix for multiple outcomes
  n.outcomes <- length(outcomes)
  if (n.outcomes > 1 && validation.results$valid) {

    if (n.outcomes == 2) {
      # 2x2 matrix: positive definite if diagonal = 1 and |off-diagonal| < 1
      # Already checked individual bounds above

    } else if (n.outcomes == 3) {
      # 3x3 matrix: create and check eigenvalues
      corr.matrix <- matrix(c(1, rho.OS.PFS, rho.OS.OR,
                              rho.OS.PFS, 1, rho.PFS.OR,
                              rho.OS.OR, rho.PFS.OR, 1),
                            nrow = 3, ncol = 3)

      eigenvalues <- eigen(corr.matrix, only.values = TRUE)$values
      validation.results$eigenvalues <- eigenvalues

      if (!all(eigenvalues > 1e-8)) {  # Small tolerance for numerical precision
        validation.results$errors <- c(validation.results$errors,
                                       sprintf("Correlation matrix is not positive definite. Eigenvalues: %s",
                                               paste(round(eigenvalues, 6), collapse = ", ")))
        validation.results$valid <- FALSE
      }
    }
  }

  # Record computation time
  end.time <- Sys.time()
  validation.results$computation.time <- as.numeric(difftime(end.time, start.time, units = "secs"))

  return(validation.results)
}
