% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorrBounds.R
\name{CorrBounds}
\alias{CorrBounds}
\title{Calculate and Validate True Frechet-Hoeffding Correlation Bounds}
\usage{
CorrBounds(
  outcomes,
  lambda.OS = NULL,
  lambda.PFS = NULL,
  p.OR = NULL,
  rho.OS.PFS = 0,
  rho.OS.OR = 0,
  rho.PFS.OR = 0
)
}
\arguments{
\item{outcomes}{Character vector of outcomes ('OS', 'PFS', 'OR')}

\item{lambda.OS}{Rate parameter for OS exponential distribution
(calculated as log(2)/median)}

\item{lambda.PFS}{Rate parameter for PFS exponential distribution
(calculated as log(2)/median)}

\item{p.OR}{Probability parameter for OR Bernoulli distribution (0 < p.OR < 1)}

\item{rho.OS.PFS}{Correlation coefficient between OS and PFS (default: 0)}

\item{rho.OS.OR}{Correlation coefficient between OS and OR (default: 0)}

\item{rho.PFS.OR}{Correlation coefficient between PFS and OR (default: 0)}
}
\value{
A list containing:
  \item{valid}{Logical indicating if all correlations are valid}
  \item{errors}{Character vector of error messages (empty if valid)}
  \item{bounds}{List of true Frechet-Hoeffding bounds for each correlation pair}
  \item{eigenvalues}{Eigenvalues of correlation matrix (if applicable)}
  \item{computation.time}{Time taken for bounds calculation}
}
\description{
This function calculates the true Frechet-Hoeffding bounds using analytical
integration for different types of variable combinations and validates
whether specified correlation coefficients are within these bounds.
}
\details{
The function calculates true Frechet-Hoeffding bounds using analytical integration:

\strong{Continuous-Continuous (OS-PFS):}
For exponential distributions, the bounds are calculated by finding the
correlation of the Frechet-Hoeffding copulas W (lower bound) and M (upper bound)
applied to exponential marginals using R's integrate() function.

\strong{Continuous-Binary (OS-OR, PFS-OR):}
The bounds are calculated using analytical integration:
\deqn{\rho_{max} = \frac{E[XY_{max}] - E[X]E[Y]}{\sigma_X \sigma_Y}}
\deqn{\rho_{min} = \frac{E[XY_{min}] - E[X]E[Y]}{\sigma_X \sigma_Y}}
where \eqn{Y_{max}} corresponds to binary=1 for largest continuous values,
and \eqn{Y_{min}} corresponds to binary=1 for smallest continuous values.

\strong{Integration Method:}
All bounds are computed using R's integrate() function for high precision
analytical integration. No grid approximation is used.
}
\examples{
# Calculate true bounds for OS, PFS, and OR
validation <- CorrBounds(
  outcomes = c('OS', 'PFS', 'OR'),
  lambda.OS = log(2)/12,  # 12-year median OS
  lambda.PFS = log(2)/6,  # 6-year median PFS
  p.OR = 0.4,             # 40\% response rate
  rho.OS.PFS = 0.6,
  rho.OS.OR = 0.3,
  rho.PFS.OR = 0.4
)

if (!validation$valid) {
  cat("Errors found:\n")
  for (error in validation$errors) {
    cat("-", error, "\n")
  }
  cat("\nTrue Frechet-Hoeffding bounds:\n")
  for (bound.name in names(validation$bounds)) {
    bound <- validation$bounds[[bound.name]]
    cat(sprintf("\%s: [\%.6f, \%.6f]\n", bound.name, bound$lower, bound$upper))
  }
}

}
