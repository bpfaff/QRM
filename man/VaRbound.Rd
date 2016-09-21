\name{VaRbound}
\alias{VaRbound}
\title{Computing lower and upper bounds for the (smallest or largest) VaR}
\description{
  \code{VaRbound()} computes lower and upper bounds for the lower or upper
  Value-at-Risk bound.
}
\usage{
VaRbound(alpha, N, qmargins, bound = c("upper", "lower"), verbose = FALSE)
}
\arguments{
  \item{alpha}{confidence level in (0,1).}
  \item{N}{tail discretization parameter; see Embrechts et al. (2013).}
  \item{qmargins}{\code{\link{list}} containing the marginal quantile
    functions.}
  \item{bound}{\code{\link{character}} string indicating the VaR bound to
    be approximated (largest (default) or smallest).}
  \item{verbose}{\code{\link{logical}} indicating whether progress
    information is displayed.}
}
\details{
  Due to the nature of the rearrangement algorithm, note that this
  purely R based implementation can be slow.
}
\value{
  \code{\link{numeric}} vector of length two, containing the lower and
  upper bound for the (chosen) Value-at-Risk estimate.
}
\author{Marius Hofert.}
\references{
  Embrechts, P., Puccetti, G., and \enc{Rüschendorf}{Rueschendorf}, L. (2013),
  Model uncertainty and VaR aggregation,
  \emph{Journal of Banking and Finance} \bold{37}(8), 2750--2764.
}
\examples{
qPar <- function(p, theta) (1-p)^(-1/theta)-1
qmar <- lapply(1:3, function(j) function(p) qPar(p, theta=2.5))
## bounds for the largest VaR
VaRbound(0.99, N=50, qmargins=qmar)
## bounds for the smallest VaR
VaRbound(0.99, N=50, qmargins=qmar, bound="lower")
}
\keyword{VaR}
\keyword{multivariate}

