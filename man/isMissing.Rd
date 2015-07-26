\name{isMissing}
\alias{isMissing}
\title{
Determine Whether Genotypes Are Missing
}
\description{
\code{isMissing} returns Boolean values indicating whether the genotypes
for a given set of samples and loci are missing from the dataset.
}
\usage{
isMissing(object, samples = Samples(object), loci = Loci(object))
}
\arguments{
  \item{object}{
An object of one of the subclasses of \code{gendata}, containing the
genotypes to be tested.
}
  \item{samples}{
A character or numeric vector indicating samples to be tested.
}
  \item{loci}{
A character or numeric vector indicating loci to be tested.
}
}
\details{
\code{isMissing} is a generic function with methods for \code{genambig}
and \code{genbinary} objects.

For each genotype in a \code{genambig} object, the function evaluates and returns
\code{Genotype(object, sample, locus)[1] == Missing(object)}.  For a
\code{genbinary} object, \code{TRUE \%in\% (Genotype(object, sample,
  locus) == Missing(object))} is returned for the genotype.  If only
one sample and locus are being evaluated, this is the Boolean value that
is returned.  If multiple samples and/or loci are being evaluated, the
function creates an array of Boolean values and recursively calls itself
to fill in the result for each element of the array.
}
\value{
  If both \code{samples} and \code{loci} are of length 1, a single
  Boolean value is returned, \code{TRUE} if the genotype is missing, and
  \code{FALSE} if it isn't.  Otherwise, the function returns a named
  array with samples in the first dimension and loci in the second
  dimension, filled with Boolean values indicating whether the genotype
  for each sample*locus combination is missing.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{Missing}}, \code{\link{Missing<-}}, \code{\link{Genotype}}, 
\code{\link{find.missing.gen}}
}
\examples{
# set up a genambig object for this example
mygen <- new("genambig", samples=c("a", "b"), loci=c("locD", "locE"))
Genotypes(mygen) <- array(list(c(122, 126), c(124, 128, 134),
                               Missing(mygen), c(156, 159)),
                          dim=c(2,2))
viewGenotypes(mygen)

# test if some individual genotypes are missing
isMissing(mygen, "a", "locD")
isMissing(mygen, "a", "locE")

# test an array of genotypes
isMissing(mygen, Samples(mygen), Loci(mygen))

}

\keyword{ manip }

