\name{Spearman}
\alias{Spearman}

\title{
  Spearman's Rank Correlation 
}

\description{
  Calculates Sperman's rank correlations. The function is a wrapper to
  \code{cor()}.  
}

\usage{
Spearman(data, ...)
}

\arguments{
  \item{data}{\code{matrix} or \code{data.frame}}
  \item{...}{ellipsis, arguments are passed down to \code{cor()}}
}

\value{
  matrix 
}

\seealso{
   \code{\link[stats]{cor}}, \code{\link{Kendall}}
}

\examples{
S <- equicorr(d = 3, rho = 0.5)
data <- rmnorm(1000, Sigma = S) 
Spearman(data) 
}

\keyword{multivariate}

