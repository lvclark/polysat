\name{simgen}
\docType{data}
\alias{simgen}
\title{Randomly Generated Data for Learning polysat}
\description{
  \code{genambig} object containing simulated data from three
  populations with 100 individuals each, at three loci.  Individuals
  are a random mixture of diploids and tetraploids.  Genotypes were
  generated according to pre-set allele frequencies.
}
\usage{data(simgen)}
\format{A \code{genambig} object with data in the \code{Genotypes},
  \code{PopInfo}, \code{PopNames}, \code{Ploidies} and \code{Usatnts}
  slots.  This is saved as an .RData file.  \code{simgen} was created
  using the code found in the \dQuote{simgen.R} file in the
  \dQuote{extdata} directory of the polysat package installation.  This
  code may be useful for inspiration on how to create a simulated dataset.}
\source{simulated data}
\seealso{\code{\link{testgenotypes}},
  \code{\linkS4class{genambig}}}
\keyword{datasets}