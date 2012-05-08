\name{QQplot}
\alias{QQplot}

\title{
  Generic Quantile-Quantile Plot 
}

\description{
  Constructs a quantile-quantile plot against a given reference
  distribution.
}

\usage{
QQplot(x, a = 0.5, reference = c("normal", "exp", "student"), ...)
}

\arguments{
  \item{x}{\code{vector}, data for QQ-plot.}
  \item{a}{\code{numeric}, the offset fraction to be used in
    \code{ppoints()}; typically in (0, 1).}
  \item{reference}{\code{character}, name of reference distribution.}
  \item{\dots}{ellipsis argument, passed down to quantile function of
  reference distribution.} 
}

\details{
  Special forms like ParetoQQ plots can also be created via this function. E.g.,
  to create a ParetoQQ plot, merely pass log(data) in place of data as the first
  parameter and use \code{reference = "exp"} as the reference
  distribution. The ParetoQQ plot should provide a linear graph when a
  log transform of the data is plotted against the exponential distribution.
}

\value{
  Produces QQ-plot and returns invisibly a list of \code{(x, y)} pairs.
}

\seealso{
  \code{\link[stats]{ppoints}}
}

\examples{
QQplot(rnorm(1000), reference = "normal") 
QQplot(rexp(1000), reference = "exp", rate = 0.3) 
}

\keyword{hplot}

