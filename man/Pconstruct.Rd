\name{Pconstruct}
\alias{Pconstruct}

\title{
  Assemble a Correlation Matrix for ML Copula Fitting 
}

\description{
  This function converts a vector of values representing the terms of a
  lower triangular matrix \eqn{A} with ones on the diagonal and
  returns the correlation matrix corresponding to the covariance matrix
  \eqn{AA'} (see page 235 in QRM).  
}

\usage{
Pconstruct(theta)
}

\arguments{
\item{theta}{\code{vector}, elements of a lower triangular matrix
  \eqn{A} with ones on the diagonal.}
}

\value{
  matrix
}

\seealso{
  \code{link{Pdeconstruct}}
}

\examples{
P <- Pconstruct(1:6) 
eigen(P) 
Pdeconstruct(P) 
}

\keyword{array}

