\name{QRM-package}
\alias{QRM-package}
\title{Quantitative Risk Modelling}
\docType{package}
\description{
  This package is designed to accompany the book Quantitative
  Risk Management: Concepts, Techniques and Tools by Alexander
  J. McNeil, Rudiger Frey and Paul Embrechts.
}

\section{Overview}{
  This package provides functions for quantitative risk management as
  introduced in the book \dQuote{Quantitative Risk Management: Concepts,
  Techniques and Tools} (henceforth: QRM). The S-Plus package
  \dQuote{QRMlib} has been made available the first author of the book
  and can be obtained by following the instructions on
  \url{http://www.ma.hw.ac.uk/~mcneil/book/QRMlib.html}. A R port of
  this package has been made available on CRAN by Scott Ulmann. However,
  the package failed the checks and hence has been moved to the
  CRAN archive (\pkg{QRMlib}, version 1.4.5.1 as of 04/25/2011). This
  package is based on \pkg{QRMlib}, but (i), not all functions have been
  ported from \pkg{QRMlib} to \pkg{QRM}, (ii) the arguments of some
  functions have been modified, and (iii) the manual pages have been
  re-ordered by topic.\cr
  A list of the not ported functions is provided in
  \code{\link{QRM-defunct}} with pointers to their replacements. This 
  was achieved by the inclusion of dependencies to the packages
  \pkg{gsl}, \pkg{numDeriv} and \pkg{timeSeries} and/or resorting to
  functions contained in the base installation of R. Second, in
  particular with respect to passing down arguments to the routines used
  in optimizations and/or argument matching, modifications to the
  functions' closures were necessary. In addition, the names of
  arguments in similar functions have been unified. Third, to provide
  the user a faster access to the manual pages of certain risk concepts,
  the functions' documentation are now ordered by concept rather than by
  the name of the functions.\cr 
  Without modifying the existing functions of \pkg{QRMlib} too much,
  neither S3- nor S4-classes and methods have been included completely
  by now in \pkg{QRM}, but the characteristic of the former package as a
  collection of functions pertinent to quantitative risk modelling have
  been kept intact. However, this might change in future releases of
  \pkg{QRM}. By now, the current package can be used almost alike
  \pkg{QRMlib}, but with the stated modifications. 
}


\references{
  McNeil, A., Frey, R. and Embrechts, P., \emph{Quantitative Risk
  Management: Concepts, Techniques and Tools}, 2005, Princeton:
  Princeton University Press.
}