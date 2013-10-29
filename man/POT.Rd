\name{POT}
\alias{POT}
\alias{fit.GPD}
\alias{findthreshold}
\alias{plotTail}
\alias{MEplot}
\alias{showRM}
\alias{xiplot}
\alias{RiskMeasures}
\alias{hill}
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
        information = c("observed", "expected"),
        optfunc = c("optim", "nlminb"), verbose = TRUE, ...)
plotTail(object, ppoints.gpd = ppoints(256), main = "Estimated tail probabilities",
         xlab = "Exceedances x", ylab = expression(1-hat(F)[n](x)), ...)
showRM(object, alpha, RM = c("VaR", "ES"),
       like.num = 64, ppoints.gpd = ppoints(256),
       xlab = "Exceedances x", ylab = expression(1-hat(F)[n](x)),
       legend.pos = "topright", ...)
findthreshold(data, ne)
MEplot(data, omit = 3., main = "Mean-Excess Plot", xlab = "Threshold",
       ylab = "Mean Excess", ...)
xiplot(data, models = 30., start = 15., end = 500., reverse = TRUE,
       ci = 0.95, auto.scale = TRUE, labels = TRUE, table = FALSE, ...)
hill(data, k, tail.index = TRUE)
hillPlot(data, option = c("alpha", "xi", "quantile"), start = 15,
         end = NA, reverse = FALSE, p = NA, ci = 0.95,
         auto.scale = TRUE, labels = TRUE, ...)
plotFittedGPDvsEmpiricalExcesses(data, threshold = NA, nextremes = NA)
RiskMeasures(out, p)
}

\arguments{
  \item{alpha}{\code{numeric}, probability level(s).}
  \item{auto.scale}{\code{logical}, whether plot should be automatically
    scaled.}
  \item{ci}{\code{numeric}, probability for asymptotic confidence
    bands.}
  \item{data}{\code{numeric}, data vector or timesSeries.}
  \item{end}{\code{integer}, maximum number of exceedances to be
    considered.}
  \item{information}{\code{character}, whether standard errors should be
    calculated with \dQuote{observed} or \dQuote{expected}
    information. This only applies to maximum likelihood type; for
    \dQuote{pwm} type \dQuote{expected} information is used if possible.}
  \item{k}{number (greater than or equal to 2) of order statistics used
           to compute the Hill plot.}
  \item{labels}{\code{logical}, whether axes shall be labelled.}
  \item{legend.pos}{if not \code{\link{NULL}}, position of \code{\link{legend}()}.}
  \item{like.num}{\code{integer}, count of evaluations of profile likelihood.}
  \item{main,xlab,ylab}{title, x axis and y axis labels.}
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
  \item{verbose}{\code{logical} indicating whether warnings are given;
    currently only applies in the case where \code{type="pwm"} and
    \code{xi > 0.5}.}
  \item{option}{\code{logical}, whether "alpha", "xi" (1 / alpha) or
    "quantile" (a quantile estimate) should be plotted.}
  \item{out}{\code{list}, returned value from fitting GPD.}
  \item{p}{\code{vector}, probability levels for risk measures.}
  \item{ppoints.gpd}{points in (0,1) for evaluating the GPD tail estimate.}
  \item{reverse}{\code{logical}, plot ordered by increasing threshold
    or number of extremes.}
  \item{RM}{\code{character}, risk measure, either "VaR" or "ES"}
  \item{start}{\code{integer}, lowest number of exceedances to be
    considered.}
  \item{table}{\code{logical}, printing of a result table.}
  \item{tail.index}{\code{logical} indicating whether the Hill estimator
    of \eqn{alpha} (the default) or \eqn{1/alpha} is computed.}
  \item{threshold}{\code{numeric}, threshold value.}
  \item{type}{\code{character}, estimation by either ML- or PWM type.}
  \item{...}{ellpsis, arguments are passed down to either \code{plot()}
    or \code{optim()} or \code{nlminb()}.}
}

\details{
  \code{hillplot()}: This plot is usually calculated from the alpha
  perspective. For a generalized Pareto analysis of heavy-tailed data
  using the gpd function, it helps to plot the Hill estimates for
  xi. See pages 286--289 in QRM. Especially note that Example 7.28
  suggests the best estimates occur when the threshold is very small,
  perhaps 0.1 of the sample size (10--50 order statistics in a sample of
  size 1000).  Hence one should NOT be using a 95 percent threshold for
  Hill estimates.\cr
  \code{MEplot()}: An upward trend in plot shows heavy-tailed
  behaviour. In particular, a straight line with positive gradient above
  some threshold is a sign of Pareto behaviour in tail. A downward trend
  shows thin-tailed behaviour whereas a line with zero gradient shows an
  exponential tail. Because upper plotting points are the average of a
  handful of extreme excesses, these may be omitted for a prettier
  plot.\cr
  \code{plotFittedGPDvsEmpiricalExcesses()}: Build a graph which plots
  the GPD fit of excesses over a threshold u and the corresponding
  empirical distribution function for observed excesses.\cr
  \code{RiskMeasures()}: Calculates risk measures (VaR or ES) based on a
  generalized Pareto model fitted to losses over a high threshold.\cr
  \code{xiplot()}: Creates a plot showing how the estimate of shape
  varies with threshold or number of extremes.
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
u <- quantile(danish, probs=0.9, names=FALSE)
plotFittedGPDvsEmpiricalExcesses(danish, threshold = u)

findthreshold(danish, 50)
mod1 <- fit.GPD(danish, threshold = u)

RiskMeasures(mod1, c(0.95, 0.99))
plotTail(mod1)

showRM(mod1, alpha = 0.99, RM = "VaR", method = "BFGS")
showRM(mod1, alpha = 0.99, RM = "ES", method = "BFGS")

mod2 <- fit.GPD(danish, threshold = u, type = "pwm")
mod3 <- fit.GPD(danish, threshold = u, optfunc = "nlminb")

## Hill plot manually constructed based on hill()

## generate data
set.seed(1)
n <- 1000 # sample size
U <- runif(n)
X1 <- 1/(1-U) # ~ F_1(x) = 1-x^{-1}, x >= 1 => Par(1)
F2 <- function(x) 1-(x*log(x))^(-1) # Par(1) with distorted SV function
X2 <- vapply(U, function(u) uniroot(function(x) 1-(x*log(x))^(-1)-u,
                                    lower=1.75, upper=1e10)$root, NA_real_)

## compute Hill estimators for various k
k <- 10:800
y1 <- hill(X1, k=k)
y2 <- hill(X2, k=k)

## Hill plot
plot(k, y1, type="l", ylim=range(y1, y2, 1),
     xlab=expression("Number"~~italic(k)~~"of upper order statistics"),
     ylab=expression("Hill estimator for"~~alpha),
     main="Hill plot") # Hill plot, good natured case (based on X1)
lines(k, y2, col="firebrick") # Hill "horror" plot (based on X2)
lines(x=c(10, 800), y=c(1, 1), col="royalblue3") # correct value alpha=1
legend("topleft", inset=0.01, lty=c(1, 1, 1), bty="n",
       col=c("black", "firebrick", "royalblue3"),
       legend=as.expression(c("Hill estimator based on"~~
                               italic(F)(x)==1-1/x,
                              "Hill estimator based on"~~
                               italic(F)(x)==1-1/(x~log~x),
                              "Correct value"~~alpha==1)))

## via hillPlot()
hillPlot(X1, option="alpha", start=10, end=800)
hillPlot(X2, option="alpha", start=10, end=800)
}

\keyword{distributions}

