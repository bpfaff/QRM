\name{Student}
\alias{Student}
\alias{dmt}
\alias{rmt}
\alias{qst}
\alias{fit.st}
\alias{fit.mst}

\title{
 Student's t Distribution 
}

\description{
  Functions for evaluating density, fitting and random variates of
  multivaraite Student's t distribution and routines for quantiles and
  fitting of univariate distribution. 
}

\usage{
dmt(x, df, mu, Sigma, log = FALSE)
rmt(n, df = 4, mu = 0, Sigma)
qst(p, mu = 0, sd = 1, df, scale = FALSE)
fit.st(data, ...)
fit.mst(data, nit = 2000, tol = 1e-10, ...)
}

\arguments{
  \item{x}{\code{matrix}, dimension \eqn{n \times d}; density is
  evaluated for each row.}
  \item{df}{\code{numeric}, degrees of freedom.}
  \item{mu}{\code{numeric}, location parameters.}
  \item{Sigma}{\code{matrix}, dispersion matrix.}
  \item{log}{\code{logical}, returning log density values.}
  \item{data}{\code{numeric}, data used for uni- and multivariate fitting.}
  \item{nit}{\code{integer}, number of iterations of EM-type algorithm.}
  \item{tol}{\code{numeric}, tolerance of improvement for stopping iteration.}
  \item{p}{\code{numeric}, probability.}
  \item{sd}{\code{numeric}, scale parameters.}
  \item{scale}{\code{logical}, scaling Student's t distribution.}
  \item{n}{\code{integer}, count of random variates.}
  \item{...}{ellipsis, arguments are passed down to \code{optim()} in
  \code{fit.st()} and to \code{MCECMupdate()} in \code{fit.mst()}.} 
}

\seealso{
  \code{link{EMupdate}}, \code{link{MCECMupdate}}, and
  \code{link{MCECM.Qfunc}}
}

\examples{
BiDensPlot(func = dmt, xpts = c(-4, 4), ypts = c(-4, 4), mu = c(0, 0),
           Sigma = equicorr(2, -0.7), df = 4)
## Quantiles of univariate Student's t
p <- c(0.90,0.95)
s <- 0.2 * 10000/sqrt(250)
qst(p, sd = s, df = 4, scale = TRUE)
## Fitting multivariate Student's t
Sigma <- diag(c(3, 4, 5)) \%*\% equicorr(3, 0.6) \%*\% diag(c(3, 4, 5)) 
mu <- c(1, 2 ,3) 
tdata <- rmt(1000, 4, mu = mu, Sigma = Sigma) 
mod1 <- fit.mst(tdata, method = "BFGS")
## DJ data
data(DJ)
r <- returns(DJ)
s <- window(r[, "MSFT"], "1993-01-01", "2000-12-31")
mod.t1 <- fit.st(100 * s)
stocks <- c("AXP","EK","BA","C","KO","MSFT",
            "HWP","INTC","JPM","DIS")
ss <- window(r[, stocks], "1993-01-01", "2000-12-31")
fridays <- time(ss)[isWeekday(time(ss), wday = 5)]
ssw <- aggregate(ss, by = fridays, FUN = sum)
mod.t2 <- fit.mst(ssw, method = "BFGS") 
}

\keyword{distribution}

