\name{ES}
\alias{ES}
\alias{ESnorm}
\alias{ESst}

\title{
  Expected Shortfall 
}

\description{
  Functions for computing the expected shortfall derived from the Normal
  or Student's t distribution (see page 45 of QRM).
}

\usage{
ESnorm(p, mu = 0, sd = 1)
ESst(p, mu = 0, sd = 1, df, scale = FALSE)
}

\arguments{
  \item{p}{\code{numeric}, probability}
  \item{mu}{\code{numeric}, location parameter}
  \item{sd}{\code{numeric}, scale parameter}
  \item{df}{\code{numeric}, degrees of freedom}
  \item{scale}{logical, scaling Student's t distribution to have variance one}
}

\value{
  numeric 
}

\examples{
p <- c(0.95, 0.99)
s <- 0.2 * 10000 / sqrt(250)
ESnorm(p)
ESst(p, sd = s, df = 4, scale = TRUE)
ESst(p, df = 4)
}

\keyword{distribution}

