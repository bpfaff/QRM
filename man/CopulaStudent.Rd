\name{CopulaStudent}
\alias{CopulaStudent}
\alias{dcopula.t}
\alias{rcopula.t}
\alias{fit.tcopula}

\title{
  Student's t Copula 
}

\description{
  Functions for copula density, generating random variates and fitting
}

\usage{
dcopula.t(Udata, df, Sigma, log = FALSE)
rcopula.t(n, df, Sigma)
fit.tcopula(Udata, method = c("all", "Kendall", "Spearman"),
            startdf = 5, ...)
}

\arguments{
  \item{df}{\code{numeric}, degrees of freedom.}
  \item{log}{\code{logical}, whether log density values should be
  returned.}
  \item{method}{\code{character}, method for fitting.}
  \item{n}{\code{integer}, count of random variates}
  \item{Sigma}{\code{matrix}, correlation matrix}
  \item{startdf}{\code{numeric}, initial DF value.}
  \item{Udata}{\code{matrix}, dimension \eqn{n \times d}, where d is the
    dimension of the copula and n is the number of pseudo-uniform
    values.}
  \item{...}{ellipsis, arguments are passed down to \code{nlminb()}.}
}

\details{
  If in the call to \code{fit.tcopula()}, \code{method = "all"}, then
  all parameters are estimated, \emph{i.e.}, the degrees of freedom and
  the dispersion parameters (initial values from Spearman
  correlations). In case of either \code{method = "Kendall"} or
  \code{method = "Spearman"}, the corresponding rank correlations are
  used and the optimization is only carried out with respect to the
  degrees of freedom parameter. The initial value for the DF is given by
  \code{startdf}. See pages 197 and 229--236 of QRM.   
}

\value{
  A vector of density values of length n for \code{dcopula.t()}. A
  matrix of random variates for \code{rcopula.t()}. A list object
  containing parameter estimates and details of fit for function
  \code{fit.tcopula()}.    
}

\seealso{
  \code{\link[stats]{nlminb}}
}

\examples{
ll <- c(0.01,0.99)
#create perspective plot for bivariate density:
BiDensPlot(func = dcopula.t, xpts = ll, ypts = ll, df = 4,
           Sigma = equicorr(2, 0.5))
S <- equicorr(d = 6, rho = 0.7)
data <- rcopula.t(2000, df = 4, Sigma = S) 
pairs(data)
## Fitting Student's Copula
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
copt1 <- fit.tcopula(Udata)
copt2 <- fit.tcopula(Udata, method = "Kendall")

}

\keyword{distribution}

