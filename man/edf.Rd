\name{edf}
\alias{edf}

\title{
  Empirical Distribution Function 
}

\description{
  This function calculates the empirical distribution function at each
  element of a vector of observations.  
}

\usage{
edf(v, adjust = FALSE)
}

\arguments{
\item{v}{\code{vector}, observations of length n.}
\item{adjust}{\code{logical}, adjustment of denominator to be (n + 1).}
}

\value{
vector 
}

\examples{
data(smi)
data(ftse100)
s1 <- window(ftse100, "1990-11-09", "2004-03-25")
s1a <- alignDailySeries(s1)
s2a <- alignDailySeries(smi)
idx <- merge(s1a, s2a)
r <-returns(idx)
rp <- series(window(r, "1994-01-01", "2003-12-31"))
rp <- rp[(rp[, 1] != 0) & (rp[, 2] !=0), ]
Udata <- apply(rp, 2, edf, adjust = 1)
plot(Udata)
}

\keyword{utilities}
