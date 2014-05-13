\name{game}
\alias{gamGPDfit}
\alias{gamGPDboot}
\title{Smooth Parameter Estimation and Bootstrapping of Generalized Pareto Distributions
  with Penalized Maximum Likelihood Estimation}
\description{
  \code{gamGPDfit()} fits the parameters of a generalized Pareto
  distribution (GPD) depending on covariates in a non- or semiparametric
  way.

  \code{gamGPDboot()} fits and bootstraps the parameters of a GPD
  distribution depending on covariates in a non- or semiparametric
  way. Applies the post-blackend bootstrap of Chavez-Demoulin and
  Davison (2005).
}
\usage{
gamGPDfit(x, threshold, nextremes = NULL, datvar, xiFrhs, nuFrhs,
          init = fit.GPD(x[,datvar], threshold = threshold,
                         type = "pwm", verbose = FALSE)$par.ests,
          niter = 32, include.updates = FALSE, eps.xi = 1e-05, eps.nu = 1e-05,
          progress = TRUE, adjust = TRUE, verbose = FALSE, \dots)
gamGPDboot(x, B, threshold, nextremes = NULL, datvar, xiFrhs, nuFrhs,
           init = fit.GPD(x[,datvar], threshold = threshold,
                          type = "pwm", verbose = FALSE)$par.ests,
           niter = 32, include.updates = FALSE, eps.xi = 1e-5, eps.nu = 1e-5,
           boot.progress = TRUE, progress = FALSE, adjust = TRUE, verbose = FALSE,
           debug = FALSE, \dots)
}
\arguments{
  \item{x}{data.frame containing the losses (in some component; can be
    specified with the argument \code{datvar}; the other components
  contain the covariates).}
  \item{B}{number of bootstrap replications.}
  \item{threshold}{threshold of the peaks-over-threshold (POT)
    method.}
  \item{nextremes}{number of excesses. This can be used to determine}
  \item{datvar}{name of the data column in \code{x} which contains the
     the data to be modeled.}
  \item{xiFrhs}{right-hand side of the formula for \eqn{\xi}{xi} in
    the \code{gam()} call for fitting \eqn{\xi}{xi}.}
  \item{nuFrhs}{right-hand side of the formula for \eqn{\nu}{nu} in
    the \code{gam()} call for fitting \eqn{\nu}{nu}.}
  \item{init}{bivariate vector containing initial values
    for \eqn{(\xi, \beta)}{(xi, beta)}.}
  \item{niter}{maximal number of iterations in the backfitting
    algorithm.}
  \item{include.updates}{\code{\link{logical}} indicating whether
  updates for xi and nu are returned as well (note: this might lead to
  objects of large size).}
  \item{eps.xi}{epsilon for stop criterion for \eqn{\xi}{xi}.}
  \item{eps.nu}{epsilon for stop criterion for \eqn{\nu}{nu}.}
  \item{boot.progress}{\code{\link{logical}} indicating
    whether progress information about \code{gamGPDboot()} is displayed.}
  \item{progress}{\code{\link{logical}} indicating whether progress information about
    \code{gamGPDfit()} is displayed. For \code{gamGPDboot()},
    \code{progress} is only passed to \code{gamGPDfit()} in the case that
    \code{boot.progress==TRUE}.}
  \item{adjust}{\code{\link{logical}} indicating whether non-real values
    of the derivatives are adjusted.}
  \item{verbose}{\code{\link{logical}} indicating whether additional
    information (in case of undesired behavior) is printed. For \code{gamGPDboot()},
    \code{progress} is only passed to \code{gamGPDfit()} if
    \code{boot.progress==TRUE}.}
  \item{debug}{\code{\link{logical}} indicating whether initial fit
  (before the bootstrap is initiated) is saved.}
  \item{\dots}{additional arguments passed to \code{gam()} (which is
    called internally; see the source code of \code{gamGPDfitUp()}).}
}
\value{
  \code{gamGPDfit()} returns either an empty list (\code{list()}; in
  case at least one of the two \code{\link{gam}()} calls in the internal
  function \code{gamGPDfitUp()} fails) or a list with the components
  \describe{
    \item{\code{xi}:}{estimated parameters \eqn{\xi}{xi};}
    \item{\code{beta}:}{estimated parameters \eqn{\beta}{beta};}
    \item{\code{nu}:}{estimated parameters \eqn{\nu}{nu};}
    \item{\code{se.xi}:}{standard error for \eqn{\xi}{xi} ((possibly
      adjusted) second-order derivative of the reparameterized
      log-likelihood with respect to \eqn{\xi}{xi}) multiplied by -1;}
    \item{\code{se.nu}:}{standard error for \eqn{\nu}{nu} ((possibly
      adjusted) second-order derivative of the reparameterized
      log-likelihood with respect to \eqn{\nu}{nu}) multiplied by -1;}
    \item{\code{xi.covar}:}{(unique) covariates for \eqn{\xi}{xi};}
    \item{\code{nu.covar}:}{(unique) covariates for \eqn{\nu}{nu};}
    \item{\code{covar}:}{available covariate combinations used for
      fitting \eqn{\beta}{beta} (\eqn{\xi}{xi}, \eqn{\nu}{nu});}
    \item{\code{y}:}{vector of excesses (exceedances minus threshold);}
    \item{\code{res}:}{residuals;}
    \item{\code{MRD}:}{mean relative distances between for all
      iterations, calculated between old parameters \eqn{(\xi, \nu)}{(xi,
	nu)} (from the last iteration) and new parameters (currently
      estimated ones);}
    \item{\code{logL}:}{log-likelihood at the estimated parameters;}
    \item{\code{xiObj}:}{\R object of type \code{gamObject} for estimated
    \eqn{\xi}{xi} (returned by \code{mgcv::gam()});}
    \item{\code{nuObj}:}{\R object of type \code{gamObject} for estimated
    \eqn{\nu}{nu} (returned by \code{mgcv::gam()});}
    \item{\code{xiUpdates}:}{if \code{include.updates} is
      \code{\link{TRUE}}, updates for \eqn{\xi}{xi} for each
      iteration. This is a list of \R objects of type \code{gamObject}
      which contains \code{xiObj} as last element;}
    \item{\code{nuUpdates}:}{if \code{include.updates} is
      \code{\link{TRUE}}, updates for \eqn{\nu}{nu} for each
      iteration. This is a list of \R objects of type \code{gamObject}
      which contains \code{nuObj} as last element;}
  }

  \code{gamGPDboot()} returns a list of length \code{B+1} where
  the first component contains the results of
  the initial fit via \code{gamGPDfit()} and the other \code{B}
  components contain the results for each replication of the
  post-blackend bootstrap. Components for which \code{\link{gam}()}
  fails (e.g., due to too few data) are given as empty lists (\code{list()}).
}
\details{
  \code{gamGPDfit()} fits the parameters \eqn{\xi}{xi} and
  \eqn{\beta}{beta} of the generalized Pareto distribution
  \eqn{\mathrm{GPD}(\xi,\beta)}{GPD(xi,beta)} depending on covariates in
  a non- or semiparametric way. The distribution function is given by
  \deqn{G_{\xi,\beta}(x)=1-(1+\xi x/\beta)^{-1/\xi},\quad
    x\ge0,}{G[xi,beta](x)=1-(1+xi x/beta)^(-1/xi), x>=0,}
  for \eqn{\xi>0}{xi>0} (which is what we assume) and
  \eqn{\beta>0}{beta>0}. Note that \eqn{\beta} is also denoted by
  \eqn{\sigma} in this package. Estimation of \eqn{\xi}{xi}
  and \eqn{\beta}{beta} by \code{gamGPDfit()} is done via penalized maximum
  likelihood estimation, where the estimators are computed with a
  backfitting algorithm. In order to guarantee convergence of this
  algorithm, a reparameterization of \eqn{\beta}{beta} in terms of the parameter
  \eqn{\nu}{nu} is done via
  \deqn{\beta=\exp(\nu)/(1+\xi).}{beta=exp(nu)/(1+xi).}
  The parameters \eqn{\xi}{xi} and \eqn{\nu}{nu} (and thus
  \eqn{\beta}{beta}) are allowed to depend on covariates (including
  time) in a non- or semiparametric way, for example:
  \deqn{\xi=\xi(\bm{x},t)=\bm{x}^{\top}\bm{\alpha}_{\xi}+h_{\xi}(t),}{xi=xi(x,t)=x^Talpha[xi]+h[xi](t),}
  \deqn{\nu=\nu(\bm{x},t)=\bm{x}^{\top}\bm{\alpha}_{\nu}+h_{\nu}(t),}{nu=nu(x,t)=x^Talpha[nu]+h[nu](t),}
  where \eqn{\bm{x}}{x} denotes the vector of covariates,
  \eqn{\bm{\alpha}_{\xi}}{alpha[xi]}, \eqn{\bm{\alpha}_{\nu}}{alpha[nu]}
  are parameter vectors and \eqn{h_{\xi}}{h[xi]}, \eqn{h_{\nu}}{h[nu]} are
  regression splines. For more details, see the references and the source
  code.

  \code{gamGPDboot()} first fits the GPD parameters via
  \code{gamGPDfit()}. It then conducts the post-blackend bootstrap of
  Chavez-Demoulin and Davison (2005). To this end, it computes the
  residuals, resamples them (\code{B} times), reconstructs the
  corresponding excesses, and refits the GPD parameters via
  \code{gamGPDfit()} again.

  Note that if \code{\link{gam}()} fails in \code{gamGPDfit()} or the
  fitting or one of the bootstrap replications in \code{gamGPDboot()},
  then the output object contains (an) empty (sub)list(s). These
  failures typically happen for too small sample sizes.
}
\author{Marius Hofert, Valerie Chavez-Demoulin.}
\references{
  Chavez-Demoulin, V., and Davison, A. C. (2005),
  Generalized additive models for sample extremes,
  \emph{Applied Statistics} \bold{54}(1), 207--222.

  Chavez-Demoulin, V., Embrechts, P., and Hofert, M.,
  An extreme value approach for modeling Operational
  Risk losses depending on covariates.
}
\examples{
## generate an example data set
years <- 2003:2012 # years
nyears <- length(years)
n <- 250 # sample size for each (different) xi
u <- 200 # threshold
rGPD <- function(n, xi, beta) ((1-runif(n))^(-xi)-1)*beta/xi # sampling GPD

set.seed(17) # setting seed
xi.true.A <- seq(0.4, 0.8, length=nyears) # true xi for group "A"
## generate losses for group "A"
lossA <- unlist(lapply(1:nyears,
                       function(y) u + rGPD(n, xi=xi.true.A[y], beta=1)))
xi.true.B <- xi.true.A^2 # true xi for group "B"
## generate losses for group "B"
lossB <- unlist(lapply(1:nyears,
                       function(y) u + rGPD(n, xi=xi.true.B[y], beta=1)))
## build data frame
time <- rep(rep(years, each=n), 2) # "2" stands for the two groups
covar <- rep(c("A","B"), each=n*nyears)
value <- c(lossA, lossB)
x <- data.frame(covar=covar, time=time, value=value)

## fit
eps <- 1e-3 # to decrease the run time for this example
require(mgcv) # due to s()
fit <- gamGPDfit(x, threshold=u, datvar="value", xiFrhs=~covar+s(time)-1,
                 nuFrhs=~covar+s(time)-1, eps.xi=eps, eps.nu=eps)
## note: choosing s(..., bs="cr") will fit cubic splines

## grab the fitted values per group and year
xi.fit <- fitted(fit$xiObj)
xi.fit. <- xi.fit[1+(0:(2*nyears-1))*n] # pick fit for each group and year
xi.fit.A <- xi.fit.[1:nyears] # fit for "A" and each year
xi.fit.B <- xi.fit.[(nyears+1):(2*nyears)] # fit for "B" and each year

## plot the fitted values of xi and the true ones we simulated from
par(mfrow=c(1,2))
plot(years, xi.true.A, type="l", ylim=range(xi.true.A, xi.fit.A),
     main="Group A", xlab="Year", ylab=expression(xi))
points(years, xi.fit.A, type="l", col="red")
legend("topleft", inset=0.04, lty=1, col=c("black", "red"),
       legend=c("true", "fitted"), bty="n")
plot(years, xi.true.B, type="l", ylim=range(xi.true.B, xi.fit.B),
     main="Group B", xlab="Year", ylab=expression(xi))
points(years, xi.fit.B, type="l", col="blue")
legend("topleft", inset=0.04, lty=1, col=c("black", "blue"),
       legend=c("true", "fitted"), bty="n")
}
\keyword{generalized Pareto distribution}
\keyword{distribution}
\keyword{multivariate}

