\name{Gauss}
\alias{Gauss}
\alias{dmnorm}
\alias{rmnorm}
\alias{fit.norm}
\alias{MardiaTest}
\alias{jointnormalTest}

\title{
  Multivariate Gauss Distribution  
}

\description{
  Functions for evaluating multivariate normal density, generating
  random variates, fitting and testing. 
}

\usage{
dmnorm(x, mu, Sigma, log = FALSE)
fit.norm(data)
rmnorm(n, mu = 0, Sigma)
MardiaTest(data)
jointnormalTest(data, dist = c("chisquare", "beta"), plot = TRUE)
}

\arguments{
  \item{data}{\code{matrix}, data set.}
  \item{dist}{\code{character}, \dQuote{chisquare} performs test against
    \eqn{\chi^2}{chi-squared} distribution, which is an approximation;
    \dQuote{beta} performs a test against a scaled beta distribution.}
  \item{log}{\code{logical}, whether log density values shall be
    returned.}
  \item{n}{\code{integer}, count of random variates.}
  \item{mu}{\code{numeric}, location parameters.}
  \item{plot}{\code{logical}, whether test result shall be plotted.} 
  \item{Sigma}{\code{matrix}, covariance matrix.}
  \item{x}{\code{matrix}, density is evaluated per row.}
}

\examples{
BiDensPlot(func = dmnorm, mu = c(0, 0), Sigma = equicorr(2, -0.7))
S <- equicorr(d = 3, rho = 0.7)
data <- rmnorm(1000, Sigma = S)
fit.norm(data)
S <- equicorr(d = 10, rho = 0.6)
data <- rmnorm(1000, Sigma = S) 
MardiaTest(data)
## Dow Jones Data
data(DJ)
r <- returns(DJ) 
stocks <- c("AXP","EK","BA","C","KO","MSFT",
            "HWP","INTC","JPM","DIS")
ss <- window(r[, stocks], "1993-01-01", "2000-12-31")
jointnormalTest(ss) 
}

\keyword{distribution}

