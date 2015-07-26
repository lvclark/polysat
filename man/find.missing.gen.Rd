\name{find.missing.gen}
\alias{find.missing.gen}
\title{
Find Missing Genotypes
}
\description{
This function returns a data frame listing the locus and sample names of
all genotypes with missing data.
}
\usage{
find.missing.gen(object, samples = Samples(object),
                 loci = Loci(object))
}
\arguments{
  \item{object}{
    A \code{genambig} or \code{genbinary} object containing the genotypes of
    interest.
}
  \item{samples}{
A character vector of all samples to be searched.  Must be a subset of
\code{Samples(object)}.
}
  \item{loci}{
A character vector of all loci to be searched.  Must be a subset of
\code{Loci(object)}.
}
}

\value{
A data frame with no row names.  The first column is named
\dQuote{Locus} and the second column is named \dQuote{Sample}.  Each
row represents one missing genotype, and gives the locus and sample
of that genotype.
}

\author{
Lindsay V. Clark
}


\seealso{
\code{\link{isMissing}}
}
\examples{
# set up the genotype data
samples <- paste("ind", 1:4, sep="")
samples
loci <- paste("loc", 1:3, sep="")
loci
testgen <- new("genambig", samples = samples, loci = loci)
Genotypes(testgen, loci="loc1") <- list(c(-9), c(102,104),
                                        c(100,106,108,110,114),
                                        c(102,104,106,110,112))
Genotypes(testgen, loci="loc2") <- list(c(77,79,83), c(79,85), c(-9),
                                        c(83,85,87,91))
Genotypes(testgen, loci="loc3") <- list(c(122,128), c(124,126,128,132),
                                        c(120,126), c(124,128,130))

# look up which samples*loci have missing genotypes
find.missing.gen(testgen)

}
\keyword{ manip }

