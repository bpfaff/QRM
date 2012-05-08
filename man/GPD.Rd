\name{GPD}
\alias{GPD}
\alias{pGPD}
\alias{qGPD}
\alias{dGPD}
\alias{rGPD}

\title{
  Generalized Pareto Distribution 
}

\description{
  Density, quantiles, and cumulative probability of the Generalized
  Pareto distribution. 
}

\usage{
pGPD(q, xi, beta = 1) 
qGPD(p, xi, beta = 1) 
dGPD(x, xi, beta = 1, log = FALSE) 
rGPD(n, xi, beta = 1)
}

\arguments{
  \item{beta}{\code{numeric}, scale parameter.}
  \item{log}{\code{logical}, whether log values of density should be
  returned.}
  \item{n}{\code{integer}, count of random variates.}
  \item{p}{\code{vector}, probabilities.}
  \item{q}{\code{vector}, quantiles.}
  \item{x}{\code{vector}, values to evaluate density.}
  \item{xi}{\code{numeric}, shape parameter.}
}

\value{
  numeric, probability (pGPD), quantile (qGPD), density (dGPD) or random
  variates (rGPD) for the GPD with scale parameter \eqn{\beta}{beta} and
  shape parameter \eqn{\xi}{xi}. 
}

\seealso{
   \code{\link{GEV}}, \code{\link{POT}}
}

\keyword{distribution}

