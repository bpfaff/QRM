\name{NH}
\alias{NH}
\alias{fit.NH}
\alias{fit.mNH}
\alias{MCECMupdate}
\alias{MCECM.Qfunc}
\alias{EMupdate}

\title{
  Normal Inverse Gaussian and Hyperbolic Distribution 
}

\description{
  Functions for fitting uni- and multivariate NIG and HYP distribution. 
}

\usage{
fit.NH(data, case = c("NIG", "HYP"), symmetric = FALSE,
       se = FALSE, ...)
fit.mNH(data, symmetric = FALSE, case = c("NIG", "HYP"),
        kvalue = NA, nit = 2000, tol = 1e-10, ...)
MCECMupdate(data, mix.pars, mu, Sigma, gamma, optpars, optfunc,
xieval=FALSE, ...)
MCECM.Qfunc(lambda, chi, psi, delta, eta, xi)
EMupdate(data, mix.pars, mu, Sigma, gamma, symmetric, 
         scaling = TRUE, kvalue = 1)
}

\arguments{
  \item{case}{\code{character}, whether NIG or HYP shall be used.}
  \item{chi}{\code{numeric}, chi parameter.}
  \item{data}{\code{numeric}, data.}
  \item{delta}{\code{numeric}, delta parameter.}
  \item{eta}{\code{numeric}, eta parameter.}
  \item{kvalue}{\code{numeric}, value to which the determinant of the
    dispersion matrix is constrained.}
  \item{lambda}{\code{numeric}, lambda parameter.}
  \item{mix.pars}{\code{vector}, values of lambda, chi and psi.}
  \item{mu}{\code{numeric}, value of location parameters.}
  \item{nit}{\code{integer}, maximum number of iterations.}
  \item{optpars}{\code{vector}, parameters to optimize over.}
  \item{optfunc}{\code{function}, the function to be optimized.}
  \item{psi}{\code{numeric}, pi parameter.}
  \item{scaling}{\code{logical}, whether determinant scaling of Sigma
    shall be fixed.}
  \item{se}{\code{logical}, whether standard errors should be
    calculated.}
  \item{Sigma}{\code{matrix}, value of Sigma.}
  \item{symmetric}{\code{logical}, whether symmetric case should be
    fitted.}
  \item{tol}{\code{numeric}, tolerance for convergence.}
  \item{gamma}{\code{numeric}, value of gamma}
  \item{xi}{\code{numeric}, xi parameter.}
  \item{xieval}{\code{logical}, whether log moment xi shall be
    evaluated.}
  \item{...}{ellipsis, arguments are passed down to \code{optim()}.}
}


\details{
  \code{fit.NH()}: See pages 78--80 of QRM. Case \sQuote{NIG} sets
  \eqn{\lambda = -1/2}{lambda = -1/2}; case \sQuote{HYP} sets
  \eqn{\lambda = 1}{lambda = 1}.\cr
  \code{fit.mNH()}: Fitting is accomplished by using a variant of the EM
  algorithm (see pages 81--83 in QRM).\cr 
  \code{MCECMupdate()}: updates estimates of mixing parameters in EM
  estimation of generalized hyperbolic (see Algorithm 3.14, steps (5)
  and (6) on page 83 in QRM).\cr
  \code{MCECM.Qfunc()}: a functional form that must be optimized when
  fitting members of generalized hyperbolic family with an MCECM
  algorithm (see function Q2 on page 82 of QRM).\cr
  \code{EMupdate()}: updates estimates of location (\eqn{\mu}{mu}),
  dispersion (\eqn{\Sigma}{Sigma}) and skewness (\eqn{\gamma}{gamma})
  parameters in EM estimation of multivariate generalized hyperbolic
  distributions (see pages 81--83 in QRM; in that case k is the
  determinant of the sample covariance matrix. \dQuote{EM} is an acronym
  for for \dQuote{Expectation-Maximization} type of algorithm
  used to fit proposed multivariate hyperbolic models to actual data). 
}

\examples{
data(DJ)
r <- returns(DJ) 
s <- window(r[, "MSFT"], "1993-01-01", "2000-12-31")
mod.NIG <- fit.NH(100 * s, method = "BFGS")
## multivariate
stocks <- c("AXP","EK","BA","C","KO","MSFT",
            "HWP","INTC","JPM","DIS")
ss <- window(r[, stocks], "1993-01-01", "2000-12-31")
fridays <- time(ss)[isWeekday(time(ss), wday = 5)]
ssw <- aggregate(ss, by = fridays, FUN = sum)
mod.mNIG <- fit.mNH(ssw, symmetric = FALSE, case = "NIG") 
}

\keyword{distribution}

