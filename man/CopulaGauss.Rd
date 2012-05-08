\name{CopulaGauss}
\alias{dcopula.gauss}
\alias{rcopula.gauss}
\alias{fit.gausscopula}

\title{
  Gauss Copula 
}

\description{
  Functions for evaluating the Gauss copula, generating random variates
  and fitting. 
}

\usage{
dcopula.gauss(Udata, Sigma, log = FALSE)
rcopula.gauss(n, Sigma)
fit.gausscopula(Udata, ...)
}

\arguments{
  \item{log}{\code{logical}, whether log density values should be
    returned.}
  \item{n}{\code{integer}, count of random variates}
  \item{Sigma}{\code{matrix}, correlation matrix.}
  \item{Udata}{\code{matrix}, pseudo-uniform data where rows are vector
    observations with all values in unit interval.}
  \item{...}{ellipsis argument, passed down to \code{nlminb()} used in
  optimization.} 
}

\value{
  For \code{dcopula.gauss()} a vector of density values of length n. For
  \code{rcopula.gauss()} a \eqn{n \times d} matrix of random variates
  and for \code{fit.gausscopula()} a list with the optimization results.    
}

\seealso{
  \code{\link[stats]{nlminb}}
}

\examples{
ll <- c(0.01,0.99)
BiDensPlot(func = dcopula.gauss, xpts = ll, ypts = ll,
           Sigma = equicorr(2, 0.5))
data <- rcopula.gauss(2000, Sigma = equicorr(d = 6, rho = 0.7)) 
pairs(data)
## Fitting Gauss Copula
data(smi)
data(ftse100)
s1 <- window(ftse100, "1990-11-09", "2004-03-25")
s1a <- alignDailySeries(s1)
s2a <- alignDailySeries(smi)
idx <- merge(s1a, s2a)
r <-returns(idx)
rp <- series(window(r, "1994-01-01", "2003-12-31"))
rp <- rp[(rp[, 1] != 0) & (rp[, 2] !=0), ]
Udata <- apply(rp, 2, edf, adjust = 1)
copgauss <- fit.gausscopula(Udata) 
}

\keyword{distribution}

