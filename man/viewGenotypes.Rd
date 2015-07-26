\name{viewGenotypes}
\alias{viewGenotypes}
\title{
Print Genotypes to the Console
}
\description{
\code{viewGenotypes} prints a tab-delimited table of samples, loci, and
alleles to the console so that genotypes can be easily viewed.
}
\usage{
viewGenotypes(object, samples = Samples(object), loci = Loci(object))
}

\arguments{
  \item{object}{
An object of
  one of the \code{gendata} subclasses, containing genotypes to be
  viewed.
}
  \item{samples}{
A numerical or character vector indicating which samples to display.
}
  \item{loci}{
A numerical or character vector indicating which loci to display.
}
}
\details{
\code{viewGenotypes} is a generic function with methods for
\code{genambig} and \code{genbinary} objects.

For a \code{genambig} object, a header line indicating sample, locus,
and allele columns is printed.
Genotypes are printed below this.  Genotypes are ordered first by locus
and second by sample.

For a \code{genbinary} object, the presence/absence matrix is printed,
organized by locus.  After the matrix for one locus is printed, a blank
line is inserted and the matrix for the next locus is printed.
}
\value{
No value is returned.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{Genotypes}}
}
\examples{
# create a dataset for this example
mygen <- new("genambig", samples=c("ind1", "ind2", "ind3", "ind4"),
             loci=c("locA", "locB"))
Genotypes(mygen) <- array(list(c(98, 104, 108), c(100, 104, 110, 114),
                               c(102, 108, 110), Missing(mygen),
                               c(132, 135), c(138, 141, 147),
                               c(135, 141, 144), c(129, 150)),
                          dim=c(4,2))

# view the genotypes
viewGenotypes(mygen)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }

