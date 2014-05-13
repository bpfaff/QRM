\name{PointProcess}
\alias{PointProcess}
\alias{extremalPP}
\alias{unmark}
\alias{fit.POT}
\alias{fit.sePP}
\alias{fit.seMPP}
\alias{stationary.sePP}
\alias{plot.MPP}
\alias{plot.PP}
\alias{plot.sePP}
\alias{sePP.negloglik}
\alias{seMPP.negloglik}
\alias{volfunction}
\alias{SEprocExciteFunc}


\title{Point Processes}

\description{
  Functions for estimating point processes.
}

\usage{
extremalPP(data, threshold = NA, nextremes = NA, ...)
unmark(PP)
fit.POT(PP, markdens = "GPD", ...)
fit.sePP(PP, model = c("Hawkes", "ETAS"), mark.influence = TRUE,
         std.errs = FALSE, ...)
fit.seMPP(PP, markdens = "GPD", model = c("Hawkes", "ETAS"),
          mark.influence = TRUE, predictable = FALSE,
          std.errs = FALSE, ...)
stationary.sePP(sePP)
sePP.negloglik(theta, PP, case)
seMPP.negloglik(theta, PP, case, markdens)
volfunction(anytimes, times, marks, theta, model)
\method{plot}{MPP}(x, ...)
\method{plot}{PP}(x, ...)
\method{plot}{sePP}(x, ...)
}

\arguments{
  \item{anytimes}{\code{vector}, times at which to calculate
    self-excitement function.}
  \item{data}{\code{timeSeries} object or \code{vector}.}
  \item{case}{\code{numeric}, indicating Hawkes or ETAS models and
  whether marks may have an influence on future points.}
  \item{markdens}{\code{character}, name of density of mark
    distribution, currently only "GPD".}
  \item{mark.influence}{\code{logical}, whether marks of marked point
    process may influence the self-excitement.}
  \item{marks}{\code{vector}, marks associated with point events.}
  \item{model}{\code{character}, name of self-exciting model.}
  \item{nextremes}{\code{integer}, count of upper extremes to be used.}  
  \item{PP}{\code{list}, a point process object of class \code{PP} or
    \code{MPP}.}
  \item{predictable}{\code{logical}, whether previous events may
    influence the scaling of mark distribution.}
  \item{sePP}{\code{list}, a fitted self-exciting process created with
  \code{fit.sePP()} or a marked self-exciting process created with
  \code{fit.seMPP()}.}
  \item{std.errs}{\code{logical}, whether standard errors should be
    computed.}
  \item{theta}{\code{vector}, parameters of self-excitement function.}
  \item{threshold}{\code{numeric}, threshold value.}
  \item{times}{\code{vector}, times of point events.}
  \item{x}{\code{list}, a (un/marked) point process object of class
    \code{PP}/\code{MPP}.}
  \item{...}{ellipsis, arguments passed to \code{plot()} or to
  \code{fit.GPD()} for \code{fit.POT()} or to \code{nlminb()} for
  functions \code{fit.sePP()} and \code{fit.seMPP} or to \code{julian()}
  for \code{extremalPP}.} 
}

\details{
  \code{extremalPP()}: returns a list describing a marked point process
  (see pages 298-301 of QRM).\cr
  \code{fit.POT()}: fits the POT (peaks-over-threshold) model to a point
  process of class \code{PP} or \code{MPP}. Note that if point process
  is of class \code{PP}, then function simply esitmates the rate of a
  homogeneous Poisson process (see pages 301--305 of QRM).\cr
  \code{fit.seMPP()}: fits a marked self-exciting process to a point
  process object of class \code{MPP}.\cr
  \code{fit.sePP()}: fits self-exciting process to a point process
  object of class \code{PP} (unmarked) or \code{MPP} (marked).\cr
  \code{seMPP.negloglik()}: evaluates negative log-likelihood of a
  marked self-exciting point process model; this objective function will
  be passed to the optimizing function.\cr
  \code{sePP.negloglik()}: evaluates negative log-likelihood of a
  self-exciting point  process model (unmarked).\cr 
  \code{stationary.sePP()}: checks a sufficient condition for
  stationarity of a self-exciting model and gives information about
  cluster size.\cr 
  \code{unmark()}: strips marks from a marked point process.\cr
  \code{volfunction()}: calculates a self-excitement function for use in
  the negloglik methods used in \code{fit.sePP()} and
  \code{fit.seMPP()}. 
}

\value{
  The function \code{extremalPP()} returns a list describing class MPP
  (marked point process) consisting of times and magnitudes of threshold
  exceedances: 
  \item{times}{vector of julian day counts (since 1/1/1960) for each
    exceedance}
  \item{marks}{vector of exceedances values (differences between value
    and threshold at each mark)}
  \item{startime}{the julian count one day prior to the first date in
    the entire timeSeries}
  \item{endtime}{value of last julian count in entire timeSeries}
  \item{threshold}{value of threshold above which exceedances are
    calculated}

  The functions \code{fit.POT()}, \code{fit.seMPP()}, and
  \code{fit.sePP()} return a list containing the fitted model.

  The \code{plot}-methods return invisibly the data for producing
  these. 
}

\seealso{
  \code{\link{GPD}}, \code{\link[stats]{nlminb}}
}

\examples{
## Extremal PP
data(sp500)
l <- -returns(sp500) 
lw <- window(l, start = "1995-12-31", end = end(l)) 
mod1 <- extremalPP(lw, ne = 100) 
mod1$marks[1:5]
mod1$threshold
mod2a <- fit.sePP(mod1, mark.influence = FALSE, std.errs = TRUE)
mod2b <- fit.seMPP(mod1, mark.influence = FALSE, std.errs = TRUE)
stationary.sePP(mod2b)
mod2c <- fit.POT(mod1, method = "BFGS")
plot(mod1)
plot(unmark(mod1))
plot(mod2a)
}

\keyword{models}

