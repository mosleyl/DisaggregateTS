% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chowlin.R
\name{chowlin}
\alias{chowlin}
\title{Function to do Chow-Lin temporal disaggregation from \insertCite{chow1971best;textual}{TSdisaggregation} and Litterman.}
\usage{
chowlin(Y, X, rho, aggMat, aggRatio, litterman = FALSE)
}
\arguments{
\item{Y}{The low-frequency response series (n_l x 1 matrix).}

\item{X}{The high-frequency indicator series (n x p matrix).}

\item{rho}{The AR(1) residual parameter (strictly between -1 and 1).}

\item{aggMat}{Aggregation matrix according to 'first', 'sum', 'average', 'last' (default is 'sum').}

\item{aggRatio}{Aggregation ratio e.g. 4 for annual-to-quarterly, 3 for quarterly-to-monthly (default is 4).}

\item{litterman}{TRUE to use litterman vcov. FALSE for Chow-Lin vcov. Default is FALSE.}
}
\value{
y	Estimated high-frequency response series (n x 1 matrix).

betaHat	Estimated coefficient vector (p x 1 matrix).

u_l	Estimated aggregate residual series (n_l x 1 matrix).
}
\description{
Used in disaggregation.R to find estimates given the optimal rho parameter.
}
\references{
\insertAllCited{}
}
\keyword{chow}
\keyword{disaggregation}
\keyword{lin}
\keyword{litterman}
\keyword{temporal}
