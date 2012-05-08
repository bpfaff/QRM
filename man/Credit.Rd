\name{Credit}
\alias{Credit}
\alias{cal.beta}
\alias{cal.claytonmix}
\alias{cal.probitnorm}
\alias{dclaytonmix}
\alias{pclaytonmix}
\alias{rclaytonmix}
\alias{dprobitnorm}
\alias{pprobitnorm}
\alias{rprobitnorm}
\alias{rlogitnorm}
\alias{rtcopulamix}
\alias{fit.binomial}
\alias{fit.binomialBeta}
\alias{fit.binomialLogitnorm}
\alias{fit.binomialProbitnorm}
\alias{momest}
\alias{rbinomial.mixture}

\title{
  Credit Risk Modelling 
}

\description{
  Functions for modelling credit risk:
  \itemize{
    \item Bernoulli mixture model with prescribed default and joint
    default probabilities
    \item Bernoulli mixture model with Clayton copula dependencies of
    default.
    \item Probitnormal Mixture of Bernoullis
    \item Beta-Binomial Distribution
    \item Logitnormal-Binomial Distribution
    \item Probitnormal-Binomial Distribution
    }
}

\usage{
cal.beta(pi1, pi2)
cal.claytonmix(pi1, pi2)
cal.probitnorm(pi1, pi2)
dclaytonmix(x, pi, theta) 
pclaytonmix(q, pi, theta) 
rclaytonmix(n, pi, theta)
rtcopulamix(n, pi, rho.asset, df)
dprobitnorm(x, mu, sigma) 
pprobitnorm(q, mu, sigma) 
rprobitnorm(n, mu, sigma)
rbinomial.mixture(n = 1000, m = 100,
                  model = c("probitnorm", "logitnorm", "beta"), ...)
rlogitnorm(n, mu, sigma)
fit.binomial(M, m)
fit.binomialBeta(M, m, startvals = c(2, 2), ses = FALSE, ...)
fit.binomialLogitnorm(M, m, startvals = c(-1, 0.5), ...)
fit.binomialProbitnorm(M, m, startvals = c(-1, 0.5), ...)
momest(data, trials, limit = 10)
}

\arguments{
\item{data}{\code{vector}, numbers of defaults in each time period.}
\item{df}{\code{numeric}, degree of freedom.}
\item{limit}{\code{intgeger}, maximum order of joint default probability
  to estimate.} 
\item{M}{\code{vector}, count of successes.}
\item{m}{\code{vector}, count of trials.}
\item{model}{\code{character}, name of mixing distribution.} 
\item{mu}{\code{numeric}, location parameter.}
\item{n}{\code{integer}, count of random variates.}
\item{pi}{\code{numeric}, default probability.}
\item{pi1}{\code{numeric}, default probability.}
\item{pi2}{\code{numeric}, joint default probability.}
\item{q}{\code{numeric}, values at which CDF should be evaluated.} 
\item{sigma}{\code{numeric}, scale parameter.}
\item{ses}{\code{logical}, whether standard errors should be returned.}
\item{startvals}{\code{numeric}, starting values.}
\item{theta}{\code{numeric}, parameter of distribution.}
\item{trials}{\code{vector}, group sizes in each time period.}
\item{x}{\code{numeric}, values at which density should be evaluated.}
\item{rho.asset}{\code{numeric}, asset correlation parameter.}
\item{...}{ellipsis, arguments are passed down to either mixing
  distribution or \code{nlminb()}.}
}

\details{
  \code{cal.beta()}: calibrates a beta mixture distribution on unit
  interval to give an exchangeable Bernoulli mixture model with
  prescribed default and joint default probabilities (see pages 354-355
  in QRM).\cr 
  \code{cal.claytonmix()}: calibrates a mixture distribution on unit
  interval to give an exchangeable Bernoulli mixture model with
  prescribed default and joint default probabilities. The mixture
  distribution is the one implied by a Clayton copula model of default
  (see page 362 in QRM).\cr 
  \code{cal.probitnorm()}: calibrates a probitnormal mixture
  distribution on unit interval to give an exchangeable Bernoulli
  mixture model with prescribed default and joint default probabilities
  (see page 354 in QRM).\cr
  \code{dclaytonmix()}, \code{pclaytonmix()}, \code{rclaytonmix()}:
  density, cumulative probability, and random generation for a mixture
  distribution on the unit interval which gives an exchangeable
  Bernoulli mixture model equivalent to a Clayton copula model (see page
  362 in QRM).\cr 
  \code{fit.binomial()}: fits binomial distribution by maximum
  likelihood.\cr 
  \code{dprobitnorm()}, \code{pprobitnorm()}, \code{rprobitnorm()}:
  density, cumulative probability and random number generation for
  distribution of random variable Q on unit interval such that the
  probit transform of Q has a normal distribution with parameters
  \eqn{\mu}{mu} and \eqn{\sigma}{sigma} (see pages 353-354 in QRM).\cr   
  \code{fit.binomialBeta()}: fit a beta-binomial distribution by maximum
  likelihood.\cr 
  \code{fit.binomialLogitnorm()}: fits a mixed binomial distribution
  where success probability has a logitnormal distribution. Lower and
  upper bounds for the input parameters M and m can be specified by
  means of the arguments \code{lower} and \code{upper}, which are passed to
  \code{nlminb()}. If convergence occurs at an endpoint of either limit,
  one need to reset lower and upper parameter estimators and run the
  function again.\cr 
  \code{fit.binomialProbitnorm()}: Fits a mixed binomial distribution
  where success probability has a probitnormal distribution. Lower and
  upper bounds for the input parameters M and m can be specified by
  means of the arguments \code{lower} and \code{upper}, which are passed to
  \code{nlminb()}. If convergence occurs at an endpoint of either limit,
  one need to reset lower and upper parameter estimators and run the
  function again.\cr
  \code{momest()}: calculates moment estimator of default probabilities
  and joint default probabilities for a homogeneous group. First
  returned value is default probability estimate; second value is
  estimate of joint default probability for two firms; and so on (see
  pages 375-376 in QRM).\cr
  \code{rbinomial.mixture()}: random variates from mixed binomial
  distribution (see pages 354-355 and pages 375-377 of QRM).\cr
  \code{rlogitnorm()}: Random number generation for distribution of
  random variable Q on unit interval such that the probit transform of Q
  has a normal distribution with parameters \eqn{\mu}{mu} and
  \eqn{\sigma}{sigma} (see pages 353-354 in QRM).\cr
  \code{rtcopulamix()}: random generation for mixing distribution on
  unit interval yielding Student's t copula model (see page 361 in QRM,
  exchangeable case of this model is considered). 
}

\seealso{
  \code{link[stats]{nlminb}}
}

\examples{
## calibrating models
pi.B <- 0.2
pi2.B <- 0.05 
probitnorm.pars <- cal.probitnorm(pi.B, pi2.B) 
probitnorm.pars 
beta.pars <- cal.beta(pi.B, pi2.B) 
beta.pars 
claytonmix.pars <- cal.claytonmix(pi.B, pi2.B) 
claytonmix.pars 
q <- (1:1000) / 1001 
q <- q[q < 0.25] 
p.probitnorm <- pprobitnorm(q, probitnorm.pars[1],
                            probitnorm.pars[2]) 
p.beta <- pbeta(q, beta.pars[1], beta.pars[2]) 
p.claytonmix <- pclaytonmix(q, claytonmix.pars[1],
                            claytonmix.pars[2]) 
scale <- range((1 - p.probitnorm), (1 - p.beta), (1 - p.claytonmix)) 
plot(q, (1 - p.probitnorm), type = "l", log = "y", xlab = "q", 
           ylab = "P(Q > q)",ylim=scale) 
lines(q, (1 - p.beta), col = 2) 
lines(q, (1 - p.claytonmix), col = 3) 
legend("topright", c("Probit-normal", "Beta", "Clayton-Mixture"), 
          lty=rep(1,3),col = (1:3))
## Clayton Mix
pi.B <- 0.0489603 
pi2.B <- 0.003126529 
claytonmix.pars <- cal.claytonmix(pi.B, pi2.B)
claytonmix.pars
q <- (1:1000) / 1001
q <- q[q < 0.25]
d.claytonmix <- dclaytonmix(q, claytonmix.pars[1], claytonmix.pars[2])
head(d.claytonmix)
## SP Data
data(spdata.raw) 
attach(spdata.raw) 
BdefaultRate <- Bdefaults / Bobligors 
## Binomial Model
mod1a <- fit.binomial(Bdefaults, Bobligors)
## Binomial Logitnorm Model
mod1b <- fit.binomialLogitnorm(Bdefaults, Bobligors) 
## Binomial Probitnorm Model
mod1c <- fit.binomialProbitnorm(Bdefaults, Bobligors)
## Binomial Beta Model
mod1d <- fit.binomialBeta(Bdefaults, Bobligors); 
## Moment estimates for default probabilities
momest(Bdefaults, Bobligors)
pi.B <- momest(Bdefaults, Bobligors)[1]
pi2.B <- momest(Bdefaults, Bobligors)[2]
## Probitnorm
probitnorm.pars <- cal.probitnorm(pi.B, pi2.B) 
q <- (1:1000)/1001
q <- q[ q < 0.25]
d.probitnorm <- dprobitnorm(q, probitnorm.pars[1], probitnorm.pars[2])
p <- c(0.90,0.95,0.975,0.99,0.995,0.999,0.9999,0.99999,0.999999)
sigma <- 0.2 * 10000 / sqrt(250)
VaR.t4 <- qst(p, df = 4, sd = sigma, scale = TRUE)
VaR.t4
detach(spdata.raw)
## Binomial Mixture Models
pi <- 0.04896 
pi2 <- 0.00321 
beta.pars <- cal.beta(pi, pi2)
probitnorm.pars <- cal.probitnorm(pi, pi2) 
n <- 1000 
m <- rep(500, n) 
mod2a <- rbinomial.mixture(n, m, "beta", shape1 = beta.pars[1],
                          shape2 = beta.pars[2]) 
mod2b <- rbinomial.mixture(n, m, "probitnorm",
                          mu = probitnorm.pars[1],
                          sigma = probitnorm.pars[2])
}

\keyword{models}

