\name{Kendall}
\alias{Kendall}

\title{
  Kendall's Rank Correlation 
}
\description{
  Calculates Kendall's rank correlations. The function is a wrapper to
  \code{cor()}.   
}

\usage{
Kendall(data, ...)
}

\arguments{
  \item{data}{\code{matrix} or \code{data.frame}.}
  \item{...}{ellipsis, arguments are passed down to \code{cor()}}
}

\value{
  matrix 
}

\seealso{
  \code{\link[stats]{cor}}, \code{\link{Spearman}}
}

\examples{
S <- equicorr(d = 3, rho = 0.5)
data <- rmnorm(1000, Sigma = S) 
Kendall(data) 
}

\keyword{multivariate}

