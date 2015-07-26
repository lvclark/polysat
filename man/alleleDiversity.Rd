\name{alleleDiversity}
\alias{alleleDiversity}
\title{
Retrieve and Count Unique Alleles
}
\description{
\code{alleleDiversity} returns the number of unique alleles and/or a
list of vectors of all unique alleles, indexed by locus and population.
}
\usage{
alleleDiversity(genobject, samples = Samples(genobject),
                loci = Loci(genobject), alleles = TRUE, counts = TRUE)
}
\arguments{
  \item{genobject}{
An object of the class \code{"genambig"}.
}
  \item{samples}{
Optional.  A character or numeric vector indicating samples to include
in the analysis.
}
  \item{loci}{
Optional.  A character or numeric vector indicating loci to include in
the analysis.
}
  \item{alleles}{
Boolean, indicating whether or not to return the alleles themselves.
}
  \item{counts}{
Boolean, indicating whether or not to return the number of unique alleles.
}
}
\value{
  Under default settings, a list is returned:
\item{alleles }{A two dimensional list.  The first dimension is indexed
  by population, with the additional element \dQuote{overall}
  representing the entire dataset.  The second dimension is indexed by
  locus.  Each element of the list is a vector, containing all unique
  alleles found for that population and locus.
  \code{Missing(genobject)} is not counted as an allele.}
\item{counts }{A matrix, indexed in the same way as \code{alleles}.
  Each element of the matrix is an integer indicating how many alleles
  were found at that population and locus.}

If the argument \code{alleles} or \code{counts} is set to \code{FALSE},
then only one of the above list elements is returned.
}
\author{
Lindsay V. Clark
}



\seealso{
\code{\link{genotypeDiversity}}
}
\examples{
# generate a dataset for this example
mygen <- new("genambig", samples=c("a","b","c","d"), loci=c("E","F"))
PopInfo(mygen) <- c(1,1,2,2)
Genotypes(mygen, loci="E") <- list(c(122,132),c(122,124,140),
                                   c(124,130,132),c(132,136))
Genotypes(mygen, loci="F") <- list(c(97,99,111),c(113,115),
                                   c(99,113),c(111,115))

# look at unique alleles
myal <- alleleDiversity(mygen)
myal$counts
myal$alleles
myal$alleles[["Pop1","E"]]
myal$alleles[["overall","F"]]
}
\keyword{ arith }

