\name{POT}
\alias{POT}
\alias{fit.GPD}
\alias{findthreshold}
\alias{plotTail}
\alias{MEplot}
\alias{showRM}
\alias{xiplot}
\alias{RiskMeasures}
\alias{hillPlot}
\alias{plotFittedGPDvsEmpiricalExcesses}

\title{
  Peaks-over-Threshold Method 
}

\description{
  Functions for fitting, analysing and risk measures according to POT/GPD
}

\usage{
fit.GPD(data, threshold = NA, nextremes = NA, type = c("ml", "pwm"),
        information = c("observed", "expected"), optfunc = c("optim", "nlminb"), ...)
findthreshold(data, ne)
plotTail(object, extend = 2, fineness = 1000, ...)
MEplot(data, omit = 3, labels=TRUE, ...)
xiplot(data, models = 30., start = 15., end = 500., reverse = TRUE, 
       ci = 0.95, auto.scale = TRUE, labels = TRUE, table = FALSE, ...)
hillPlot(data, option = c("alpha", "xi", "quantile"), start = 15, 
         end = NA, reverse = FALSE, p = NA, ci = 0.95,
         auto.scale = TRUE, labels = TRUE, ...)
plotFittedGPDvsEmpiricalExcesses(data, threshold = NA, nextremes = NA)
showRM(object, alpha, RM = c("VaR", "ES"), extend = 2, ci.p = 0.95,
       like.num = 50., ...)
RiskMeasures(out, p)
}

\arguments{
  \item{alpha}{\code{numeric}, probability level(s).}
  \item{auto.scale}{\code{logical}, whether plot should be automatically
    scaled.}
  \item{ci}{\code{numeric}, probability for asymptotic confidence
    bands.}
  \item{ci.p}{\code{numeric}, confidence levels.}
  \item{data}{\code{numeric}, data vector or timesSeries.}
  \item{end}{\code{integer}, maximum number of exceedances to be
    considered.}
  \item{extend}{\code{numeric}, extension of plotting range.}
  \item{fineness}{\code{integer}, count of points at which to evaluate
    the tail estimate.}
  \item{information}{\code{character}, whether standard errors should be
    calculated with \dQuote{observed} or \dQuote{expected}
    information. This only applies to maximum likelihood type; for
  \dQuote{pwm} type \dQuote{expected} information is used if possible.}
  \item{labels}{\code{logical}, whether axes shall be labelled.}
  \item{like.num}{\code{integer}, count of evaluations of profile likelihood.}
  \item{type}{\code{character}, estimation by either ML- or PWM type.}
  \item{models}{\code{integer}, count of consecutive gpd models to be
    fitted; \emph{i.e.}, the count of different thresholds at which to
    re-estimate \eqn{\xi}{xi}; this many \eqn{\xi}{xi} estimates will be
    plotted.}
  \item{ne}{\code{integer}, count of excesses above the threshold.}
  \item{nextremes}{\code{integer}, count of upper extremes to be used.}
  \item{object}{\code{list}, returned value from fitting GPD}
  \item{omit}{\code{integer}, count of upper plotting points to be
    omitted.}
  \item{optfunc}{\code{character}, function used for ML-optimization.}
  \item{option}{\code{logical}, whether "alpha", "xi" (1 / alpha) or
    "quantile" (a quantile estimate) should be plotted.}
  \item{out}{\code{list}, returned value from fitting GPD.}
  \item{p}{\code{vector}, probability levels for risk measures.}
  \item{reverse}{\code{logical}, plot ordered by increasing threshold
    or number of extremes.}
  \item{RM}{\code{character}, risk measure, either "VaR" or "ES"}
  \item{start}{\code{integer}, lowest number of exceedances to be
    considered.}
  \item{table}{\code{logical}, printing of a result table.}
  \item{threshold}{\code{numeric}, threshold value.}
  \item{...}{ellpsis, arguments are passed down to either \code{plot()}
    or \code{optim()} or \code{nlminb()}.}
}

\details{
  \code{MEplot()}: An upward trend in plot shows heavy-tailed
  behaviour. In particular, a straight line with positive gradient above
  some threshold is a sign of Pareto behaviour in tail. A downward trend
  shows thin-tailed behaviour whereas a line with zero gradient shows an
  exponential tail. Because upper plotting points are the average of a
  handful of extreme excesses, these may be omitted for a prettier
  plot.\cr
  \code{RiskMeasures()}: calculates risk measures (VaR or ES) based on a
  generalized Pareto model fitted to losses over a high threshold.\cr 
  \code{xiplot()}: creates a plot showing how the estimate of shape
  varies with threshold or number of extremes.\cr
  \code{hillplot()}: This plot is usually calculated from the alpha
  perspective. For a generalized Pareto analysis of heavy-tailed data
  using the gpd function, it helps to plot the Hill estimates for
  xi. See pages 286--289 in QRM. Especially note that Example 7.28
  suggests the best estimates occur when the threshold is very small,
  perhaps 0.1 of the sample size (10--50 order statistics in a sample of
  size 1000).  Hence one should NOT be using a 95 percent threshold for
  Hill estimates.\cr 
  \code{plotFittedGPDvsEmpiricalExcesses()}: Build a graph which plots
  the GPD fit of excesses over a threshold u and the corresponding
  empirical distribution function for observed excesses. 
}

\seealso{
  \code{\link{GEV}}
}

\examples{
data(danish)
plot(danish)
MEplot(danish)
xiplot(danish) 
hillPlot(danish, option = "alpha", start = 5, end = 250, p = 0.99)
hillPlot(danish, option = "alpha", start = 5, end = 60, p = 0.99)
plotFittedGPDvsEmpiricalExcesses(danish, nextremes = 109)
plotFittedGPDvsEmpiricalExcesses(danish, threshold = 10)
findthreshold(danish, 50)
mod1 <- fit.GPD(danish, threshold = 10) 
RiskMeasures(mod1, c(0.95, 0.99)) 
plotTail(mod1)
showRM(mod1, 0.99, RM = "VaR", method = "BFGS")
showRM(mod1, 0.99, RM = "ES", method = "BFGS")
mod2 <- fit.GPD(danish, threshold = 10, type = "pwm")
mod3 <- fit.GPD(danish, threshold = 10, optfunc = "nlminb")
}

\keyword{distributions}

