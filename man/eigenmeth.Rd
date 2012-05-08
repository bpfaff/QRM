\name{eigenmeth}
\alias{eigenmeth}

\title{
  Make Matrix Positive Definite 
}

\description{
  The function adjusts a negative definite symmetric matrix to make it
  positive definite.  
}

\usage{
eigenmeth(mat, delta = 0.001)
}

\arguments{
\item{mat}{\code{matrix}, a symmetric matrix}
\item{delta}{\code{numeric}, new size of smallest eigenvalues}
}

\value{
a positive-definite matrix 
}

\details{
See page 231 of QRM. 
}

\keyword{array}

