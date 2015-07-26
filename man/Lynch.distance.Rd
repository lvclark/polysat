\name{Lynch.distance}
\alias{Lynch.distance}
\title{
Calculate Band-Sharing Dissimilarity Between Genotypes
}
\description{
Given two genotypes in the form of vectors of unique alleles, a
dissimilarity is calculated as:
1 - (number of alleles in common)/(average number of alleles per genotype).
}
\usage{
Lynch.distance(genotype1, genotype2, usatnt = NA, missing = -9)
}
\arguments{
  \item{genotype1}{
A vector containing all alleles for a particular sample and locus.
Each allele is only present once in the vector.
}
  \item{genotype2}{
A vector of the same form as \code{genotype1}, for another sample at the
same locus.
}
\item{usatnt}{
The microsatellite repeat length for this locus (ignored by the function).
}
  \item{missing}{
The symbol used to indicate missing data in either genotype vector.
}
}
\details{
Lynch (1990) defines a simple measure of similarity between DNA
fingerprints.  This is 2 times the number of bands that two fingerprints
have in common, divided by the total number of bands that the two genotypes
have.  \code{Lynch.distance} returns a dissimilarity, which is 1 minus
the similarity.
}
\value{
If the first element of either or both input genotypes is equal to
\code{missing}, NA is returned.

Otherwise, a numerical value is returned.  This is one minus the
similarity.  The similarity is calculated as the number of alleles that
the two genotypes have in common divided by the mean length of the two
genotypes.
}
\references{
Lynch, M.  (1990)  The similarity index and DNA fingerprinting.
\emph{Molecular Biology and Evolution} \bold{7}, 478-484.
}
\author{
Lindsay V. Clark
}

\seealso{
  \code{\link{Bruvo.distance}},
  \code{\link{meandistance.matrix}}
}
\examples{
Lynch.distance(c(100,102,104), c(100,104,108))
Lynch.distance(-9, c(102,104,110))
Lynch.distance(c(100), c(100,104,106))
}
\keyword{ arith }

