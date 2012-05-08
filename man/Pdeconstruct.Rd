\name{Pdeconstruct}
\alias{Pdeconstruct}

\title{
  Disassemble a Correlation Matrix for ML Copula Fitting 
}

\description{
  This function takes a correlation matrix \eqn{P} and returns the
  elements of a lower-triangular matrix \eqn{A} with ones on the
  diagonal such that \eqn{P} is the corelation matrix corresponding to
  the covariance matrix \eqn{AA'} (see page 235 in QRM).   
}

\usage{
Pdeconstruct(P)
}

\arguments{
\item{P}{\code{matrix}, a correlation matrix}
}

\value{
  vector
}

\seealso{
  \code{\link{Pconstruct}}
}

\examples{
P <- Pconstruct(1:6) 
Pdeconstruct(P) 
}

\keyword{array}
