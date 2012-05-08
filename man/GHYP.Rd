\name{GHYP}
\alias{dghyp}
\alias{dmghyp}
\alias{dsmghyp}
\alias{dghypB}
\alias{rghyp}
\alias{rmghyp}
\alias{rghypB}

\title{
  Uni- and Multivariate Generalized Hyperbolic Distribution 
}

\description{
  Values of density and random number generation for uni- and
  multivariate Generalized Hyperbolic distribution in new QRM
  parameterization \eqn{(\chi, \psi, \gamma)}{(chi, psi, gamma)} and in
  standard parametrization \eqn{(\alpha, \beta, \delta)}{(alpha, beta,
  delta)}; univariate only. See pp. 77--81 in QRM. The special case of a
  multivariate symmetric GHYP is implemented seperately as function
  \code{dsmghyp()}.    
}

\usage{
dghyp(x, lambda, chi, psi, mu = 0, gamma = 0, log = FALSE) 
dmghyp(x, lambda, chi, psi, mu, Sigma, gamma, log = FALSE)
dsmghyp(x, lambda, chi, psi, mu, Sigma, log = FALSE)
dghypB(x, lambda, delta, alpha, beta = 0, mu = 0, log = FALSE) 
rghyp(n, lambda, chi, psi, mu = 0, gamma = 0)
rmghyp(n, lambda, chi, psi, Sigma, mu, gamma)
rghypB(n, lambda, delta, alpha, beta = 0, mu = 0)
}

\arguments{
  \item{alpha}{\code{numeric}, parameter(s).}
  \item{beta}{\code{numeric}, skewness parameter.}
  \item{chi}{\code{numeric},  mixing parameter(s).}
  \item{delta}{\code{numeric}, parameter(s).}
  \item{gamma}{\code{numeric}, skewness parameter(s).}
  \item{lambda}{\code{numeric}, mixing parameter(s).}
  \item{log}{\code{logical}, should log density be returned; default is
    \code{FALSE}.}
  \item{mu}{\code{numeric}, location parameter(s).}
  \item{n}{\code{integer}, count of random variates.}
  \item{psi}{\code{numeric}, mixing parameter(s).}
  \item{Sigma}{\code{matrix}, dispersion matrix for multivaraiate GHYP.}
  \item{x}{\code{vector}, values to evaluate density.}
}

\details{
  The univariate QRM parameterization is defined in terms of parameters
\eqn{\chi, \psi, \gamma}{chi, psi, gamma} instead of the
\eqn{\alpha, \beta, \delta}{alpha, beta, delta} model used by Blaesild
(1981). If \eqn{\gamma = 0}{gamma = 0}, a normal variance mixture where
  the mixing variable \eqn{W} has a Generalized Inverse Gaussian
  distribution (GIG) with parameters \eqn{\lambda, \chi, \psi}{lambda,
  chi, psi} is given, with heavier tails. If \eqn{\gamma > 0}{gamma >
  0}, a normal mean-variance mixture where the mean is also perturbed to
  equal \eqn{\mu + (W * \gamma)}{mu + (W * gamma)} which introduces
  asymmetry as well, is obtained. Values for \eqn{\lambda}{lambda} and
  \eqn{\mu}{mu} are identical in both QRM and B parameterizations. The
  dispersion matrix \eqn{\Sigma}{Sigma} does not appear as argument in
  the univariate case since its value is identically one.  
}

\note{
  Density values from dgyhp() should be identical to those from dghypB()
  if the \eqn{\alpha, \beta, \delta}{alpha, beta, delta} parameters of
  the B type are translated to the corresponding \eqn{\gamma, \chi,
  \psi}{gamma,chi, psi} parameters of the QRM type by formulas on pp
79--80 in QRM.\cr
If \eqn{\gamma}{gamma} is a vector of zeros, the distribution is
  elliptical and \code{dsmghyp()} is utilised in \code{dmghyp()}. If
  \eqn{\lambda = (d + 1) / 2}{lambda = (d + 1) / 2}, a d-dimensional
  hyperbolic density results. If \eqn{\lambda = 1}{lambda = 1}, the
  univariate marginals are one-dimensional hyperbolics. If \eqn{\lambda
  = -1/2}{lambda = -1/2}, the distribution is Normal Inverse Gaussian
  (NIG). If \eqn{\lambda > 0}{lambda > 0} and \eqn{\chi = 0}{chi = 0},
  one obtains a Variance Gamma distribution (VG). If one can define a
  constant \eqn{\nu}{nu} such that \eqn{\lambda = (-1/2) * \nu}{lambda =
  (-1/2) * nu} and \eqn{\chi = \nu}{chi = nu} then one obtains a
  multivariate skewed-t distribution. See p. 80 of QRM for details.  
}

\value{
  numeric, value(s) of density or log-density (dghyp, dmghyp, dsmghyp
  and dghypB) or random sample (rghyp, rmghyp, rghypB) 
}

\examples{
old.par <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
ll <- c(-4, 4)
BiDensPlot(func = dmghyp, xpts = ll, ypts = ll, mu = c(0, 0),
           Sigma = equicorr(2, -0.7), lambda = 1, chi = 1, psi = 1,
           gamma = c(0, 0))
BiDensPlot(func = dmghyp, type = "contour", xpts = ll, ypts = ll,
           mu = c(0, 0), Sigma = equicorr(2, -0.7), lambda = 1,
           chi = 1, psi = 1, gamma = c(0, 0))
BiDensPlot(func = dmghyp, xpts = ll, ypts = ll, mu = c(0, 0),
           Sigma = equicorr(2, -0.7), lambda = 1, chi = 1, psi = 1,
           gamma = c(0.5, -0.5))
BiDensPlot(func = dmghyp, type = "contour", xpts = ll, ypts = ll,
           mu = c(0, 0), Sigma = equicorr(2, -0.7), lambda = 1,
           chi = 1, psi = 1, gamma = c(0.5, -0.5))
par(old.par)
}


\keyword{distribution}

