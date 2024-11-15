% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TempDisaggDGP.R
\name{TempDisaggDGP}
\alias{TempDisaggDGP}
\title{High and Low-Frequency Data Generating Processes}
\usage{
TempDisaggDGP(
  n_l,
  n,
  aggRatio = 4,
  p = 1,
  beta = 1,
  sparsity = 1,
  method = "Chow-Lin",
  aggMat = "sum",
  rho = 0,
  mean_X = 0,
  sd_X = 1,
  sd_e = 1,
  simul = FALSE,
  sparse_option = "random",
  setSeed = 42
)
}
\arguments{
\item{n_l}{Integer. Size of the low-frequency series.}

\item{n}{Integer. Size of the high-frequency series.}

\item{aggRatio}{Integer. Aggregation ratio between low and high frequency (default is 4).}

\item{p}{Integer. Number of high-frequency indicator series to include.}

\item{beta}{Numeric. Value for the positive and negative elements of the coefficient vector.}

\item{sparsity}{Numeric. Sparsity percentage of the coefficient vector (value between 0 and 1).}

\item{method}{Character. The DGP of residuals to use ('Denton', 'Denton-Cholette', 'Chow-Lin', 'Fernandez', 'Litterman').}

\item{aggMat}{Character. Aggregation matrix type ('first', 'sum', 'average', 'last').}

\item{rho}{Numeric. Residual autocorrelation coefficient (default is 0).}

\item{mean_X}{Numeric. Mean of the design matrix (default is 0).}

\item{sd_X}{Numeric. Standard deviation of the design matrix (default is 1).}

\item{sd_e}{Numeric. Standard deviation of the errors (default is 1).}

\item{simul}{Logical. If \code{TRUE}, the design matrix and the coefficient vector are fixed (default is \code{FALSE}).}

\item{sparse_option}{Character or Integer. Option to specify sparsity in the coefficient vector ('random' or integer value). Default is "random".}

\item{setSeed}{Integer. Seed value for reproducibility when \code{simul} is set to \code{TRUE}.}
}
\value{
A list containing the following components:
\itemize{
\item \code{y_Gen}: Generated high-frequency response series (an \eqn{n \times 1} matrix).
\item \code{Y_Gen}: Generated low-frequency response series (an \eqn{n_l \times 1} matrix).
\item \code{X_Gen}: Generated high-frequency indicator series (an \eqn{n \times p} matrix).
\item \code{Beta_Gen}: Generated coefficient vector (a \eqn{p \times 1} matrix).
\item \code{e_Gen}: Generated high-frequency residual series (an \eqn{n \times 1} matrix).
}
}
\description{
This function generates a high-frequency response vector \eqn{y}, following the relationship \eqn{y = X\beta + \epsilon}, where \eqn{X} is a matrix of indicator series and \eqn{\beta} is a potentially sparse coefficient vector. The low-frequency vector \eqn{Y} is generated by aggregating \eqn{y} according to a specified aggregation method.
}
\details{
The aggregation ratio (\code{aggRatio}) determines the ratio between the low and high-frequency series (e.g., \code{aggRatio = 4} for annual-to-quarterly). If the number of observations \eqn{n} exceeds \eqn{aggRatio \times n_l}, the aggregation matrix will include zero columns for the extrapolated values.

The function supports several data generating processes (DGP) for the residuals, including 'Denton', 'Denton-Cholette', 'Chow-Lin', 'Fernandez', and 'Litterman'. These methods differ in how they generate the high-frequency data and residuals, with optional autocorrelation specified by \code{rho}.
}
\examples{
data <- TempDisaggDGP(n_l=25,n=100,p=10,rho=0.5)
X <- data$X_Gen
Y <- data$Y_Gen
}
\keyword{DGP}
\keyword{high-frequency}
\keyword{low-frequency}
\keyword{sparse}
