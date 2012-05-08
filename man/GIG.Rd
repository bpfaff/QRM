\name{GIG}
\alias{GIG}
\alias{EGIG}
\alias{ElogGIG}
\alias{rGIG}

\title{Generalized Inverse Gaussian Distribution}

\description{
  Calculates (log) moments of univariate generalized inverse Gaussian
  (GIG) distribution and generating random variates. 
}

\usage{
EGIG(lambda, chi, psi, k = 1)
ElogGIG(lambda, chi, psi)
rGIG(n, lambda, chi, psi, envplot = FALSE, messages = FALSE)
}

\arguments{
  \item{chi}{\code{numeric}, chi parameter.}
  \item{envplot}{\code{logical}, whether plot of rejection envelope
    should be created.}
  \item{k}{\code{integer}, order of moments.}
  \item{lambda}{\code{numeric}, lambda parameter.}
  \item{messages}{\code{logical}, whether a message about rejection rate
    should be returned.}
  \item{n}{\code{integer}, count of random variates.}
  \item{psi}{\code{numeric}, psi parameter.}
}

\details{
  Normal variance mixtures are frequently obtained by perturbing the
  variance component of a normal distribution; here this is done by
  multiplying the square root of a mixing variable assumed to have a GIG
  distribution depending upon three parameters \eqn{(\lambda, \chi,
    \psi)}{(lambda, chi, psi)}. See p.77 in QRM.\cr
  Normal mean-variance mixtures are created from normal variance
  mixtures by applying another perturbation of the same mixing variable
  to the mean component of a normal distribution. These perturbations
  create Generalized Hyperbolic Distributions. See pp. 78--81 in QRM. A
  description of the GIG is given on page 497 in QRM Book.
}

\value{
  (log) mean of distribution or vector random variates in case of
  \code{rgig()}.} 

\keyword{distribution}

