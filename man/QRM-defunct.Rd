\name{QRM-defunct}
\alias{QRM-defunct}
\alias{mk.returns}
\alias{plotMultiTS}
\alias{signalSeries}
\alias{aggregateQuarterlySeries}
\alias{aggregateMonthlySeries}
\alias{aggregateWeeklySeries}
\alias{aggregateSignalSeries}
\alias{ConvertDFToTimeSeries}
\alias{CovToCor}
\alias{symmetrize}
\alias{hessb}
\alias{fit.Archcopula2d}
\alias{besselM3}
\alias{psifunc}
\alias{kurtosisSPlus}
\alias{fit.tcopula.rank}
\alias{fit.GPDb}
\alias{lbeta}

\docType{data}

\title{Defunct Functions in Package QRM}

\description{
  The functions listed below which were contained in the package QRMlib
  are now defunct. The user is referred to the suggested functions as an
  alternative. 
}

\details{
  \code{aggregateMonthlySeries()} is defunct. use \code{aggregate()} in
  package \pkg{timeSeries}.\cr
  \code{aggregateQuarterlySeries} is defunct. use \code{aggregate()} in
  package \code{timeSeries}.\cr 
  \code{aggregateSignalSeries()} is defunct. use \code{aggregate()} in
  package \pkg{timeSeries}.\cr
  \code{aggregateWeeklySeries()} is defunct. use \code{aggregate()} in
  package \pkg{timeSeries}.\cr
  \code{besselM3()} is defunct. use \code{besselK()} in package
  \pkg{base}.\cr 
  \code{ConvertDFToTimeSeries()} is defunct. use \code{timeSeries()} in
  package \pkg{timeSeries}.\cr
  \code{CovToCor()} is defunct. use \code{cov2cor()} in package
  \pkg{stats}.\cr 
  \code{fit.Archcopula2d()} is defunct. use \code{fit.AC()}.\cr
  \code{fit.GPDb()} is defunct. use \code{fit.GPD()}.\cr
  \code{fit.tcopula.rank()} is defunct. use \code{fit.tcopula()}.\cr
  \code{hessb()} is defunct. use \code{hessian()} in package
  \pkg{numDeriv}.\cr 
  \code{kurtosisSPlus()} is defunct. use \code{kurtosis()} in package
  \pkg{timeDate}.\cr 
  \code{lbeta()} is defunct. use \code{lbeta()} in package
  \pkg{base}.\cr  
  \code{mk.returns()} is defunct. use \code{returnSeries()} in package
  \pkg{timeSeries}.\cr  
  \code{plotMultiTS()} is defunct. use \code{plot()} in package
  \pkg{timeSeries}.\cr
  \code{psifunc()} is defunct. use \code{psi()} in package \pkg{gsl}.\cr 
  \code{signalSeries()} is defunct. use \code{series()} in package
  \pkg{timeSeries}.\cr  
  \code{symmetrize()} is defunct. use \code{forceSymmetric()} in package
  \pkg{Matrix}.\cr 
}

\keyword{data}