\name{BiDensPlot}
\alias{BiDensPlot}

\title{
  Bivariate Density Plot 
}

\description{
  Generates eiether a perspective or a contour plot of a bivariate
  density.  
}

\usage{
BiDensPlot(func, xpts = c(-2, 2), ypts = c(-2, 2), npts = 50,
           type = c("persp", "contour"), ...)
}

\arguments{
  \item{func}{\code{function}, the name of a bivariate density
    function.}
  \item{xpts}{\code{vector}, interval of x.}
  \item{ypts}{\code{vector}, interval of y.}
  \item{npts}{\code{integer}, number of subdivision points between x and
    y over the specified range xpts to ypts.}
  \item{type}{\code{character}, the plot type, either a perspective
   or a contour plot.}
  \item{...}{\code{ellipsis}, arguments are passed to the call of
  \code{func}.}  
}

\value{
  Returns invisibly a list of \code{(x, y, z)} triplet.
}

\examples{
BiDensPlot(func = dmnorm, mu = c(0, 0), Sigma = equicorr(2, -0.7)) 
}

\keyword{hplot}
