\name{Gumbel}
\alias{dGumbel}
\alias{qGumbel}
\alias{pGumbel}
\alias{rGumbel}

\title{
  Gumbel Distribution 
}

\description{
  Density, quantiles, and cumulative probability of the Gumbel
  distribution. The standard Gumbel has \eqn{\mu}{mu} value of 0 and
  \eqn{\sigma}{sigma} value of one. 
}

\usage{
dGumbel(x, mu = 0, sigma = 1, log = FALSE) 
qGumbel(p, mu = 0, sigma = 1) 
pGumbel(q, mu = 0, sigma = 1)
rGumbel(n, mu = 0, sigma = 1)
}

\arguments{
  \item{log}{\code{logical}, whether log values of density should be
    returned.}
  \item{mu}{\code{numeric}, location parameter.}
  \item{n}{\code{integer}, count of random variates.}
  \item{p}{\code{vector}, probabilities.}
  \item{q}{\code{vector}, quantiles.}
  \item{sigma}{\code{numeric}, scale parameter.}
  \item{x}{\code{vector}, values to evaluate density.}
}

\value{
  numeric, probability (\code{pGumbel()}), quantile (\code{qGumbel()}),
  density (\code{dGumbel()}) or random variates (\code{rGumbel()}) for
  the Gumbel distribution with location parameter \eqn{\mu}{mu} and
  scale parameter \eqn{\sigma}{sigma}.  
}

\examples{
rGumbelSim <- rGumbel(1000, 1.0, 2.5)
quantValue <- 4.5
pGEV(q = quantValue, xi = 0, mu = 1.0, sigma = 2.5) 
pGumbel(q = quantValue, mu = 1.0, sigma = 2.5)
}

\keyword{distribution}
