\name{GEV}
\alias{GEV}
\alias{pGEV}
\alias{qGEV}
\alias{dGEV}
\alias{rGEV}
\alias{fit.GEV}

\title{
  Generalized Extreme Value Distribution 
}

\description{
  Density, quantiles, cumulative probability, and fitting of the Generalized
  Extreme Value distribution. 
}

\usage{
pGEV(q, xi, mu = 0, sigma = 1) 
qGEV(p, xi, mu = 0, sigma = 1) 
dGEV(x, xi, mu = 0, sigma = 1, log = FALSE) 
rGEV(n, xi, mu = 0, sigma = 1)
fit.GEV(maxima, ...)
}

\arguments{
  \item{log}{\code{logical}, whether log values of density should be returned,
    default is \code{FALSE}.}
  \item{maxima}{\code{vector}, block maxima data}
  \item{mu}{\code{numeric}, location parameter.}
  \item{n}{\code{integer}, count of random variates.}
  \item{p}{\code{vector}, probabilities.}
  \item{q}{\code{vector}, quantiles.}
  \item{sigma}{\code{numeric}, scale parameter.}
  \item{x}{\code{vector}, values to evaluate density.}
  \item{xi}{\code{numeric}, shape parameter.}
  \item{...}{ellipsis, arguments are passed down to \code{optim()}.}
}

\value{
  numeric, probability (pGEV), quantile (qGEV), density (dGEV) or random
  variates (rGEV) for the GEV distribution with shape parameter
  \eqn{\xi}{xi}, location parameter \eqn{\mu}{mu} and scale parameter
  \eqn{\sigma}{sigma}. A list object in case of \code{fit.GEV()}.
}

\seealso{
   \code{\link{GPD}}
}

\examples{
quantValue <- 4.5
pGEV(q = quantValue, xi = 0, mu = 1.0, sigma = 2.5) 
pGumbel(q = quantValue, mu = 1.0, sigma = 2.5)
## Fitting to monthly block-maxima
data(nasdaq)
l <- -returns(nasdaq)
em <- timeLastDayInMonth(time(l))
monmax <- aggregate(l, by = em, FUN = max) 
mod1 <- fit.GEV(monmax) 
}

\keyword{distribution}

