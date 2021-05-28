\name{write.ATetra}
\alias{write.ATetra}
\title{
Write Genotypes in ATetra Format
}
\description{
 \code{write.ATetra} uses genotype and population information contained
 in a \code{genambig} object to create a
text file of genotypes in the ATetra format.
}
\usage{
write.ATetra(object, samples = Samples(object),
             loci = Loci(object), file = "")
}
\arguments{
  \item{object}{
A \code{genambig} object containing the dataset of interest.  Genotypes,
population identities, population names, and the dataset description are
used for creating the file.  Ploidies must be set to 4.
}
  \item{samples}{
    A character vector of samples to write to the file.  This is a subset
    of \code{Samples(object)}.
}
  \item{loci}{
    A character vector of loci to write to the file.  This is a subset of
    \code{Loci(object)}.
}
  \item{file}{
A character string indicating the path and name to which to write the file.
}
}
\details{
Note that missing data are not allowed in ATetra, although
\code{write.ATetra} will still process missing data.  When it does so, it
leaves all alleles blank in the file for that particular sample and
locus, and also prints a warning indicating which sample and locus had
missing data.
}
\value{
A file is written but no value is returned.
}
\references{
  \url{http://www.vub.ac.be/APNA/ATetra_Manual-1-1.pdf}

  van Puyvelde, K., van Geert, A. and Triest, L.  (2010)  ATETRA, a new
  software program to analyze tetraploid microsatellite data: comparison
  with TETRA and TETRASAT.  \emph{Molecular Ecology Resources}
  \bold{10}, 331-334.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.ATetra}}, \code{\link{write.Tetrasat}},
\code{\link{write.GeneMapper}}, \code{\link{write.POPDIST}}
}
\examples{
# set up sample data (usually done by reading files)
mysamples <- c("ind1", "ind2", "ind3", "ind4")
myloci <- c("loc1", "loc2")
mygendata <- new("genambig", samples=mysamples, loci=myloci)
mygendata <- reformatPloidies(mygendata, output="one")
Genotypes(mygendata, loci="loc1") <- list(c(202,204), c(204),
                                          c(200,206,208,212),
                                          c(198,204,208))
Genotypes(mygendata, loci="loc2") <- list(c(78,81,84),
                                          c(75,90,93,96,99),
                                          c(87), c(-9))
PopInfo(mygendata) <- c(1,2,1,2)
PopNames(mygendata) <- c("this pop", "that pop")
Ploidies(mygendata) <- 4
Description(mygendata) <- "Example for write.ATetra."

\dontrun{
# write an ATetra file
write.ATetra(mygendata, file="atetratest.txt")

# view the file
cat(readLines("atetratest.txt"),sep="\n")
}
\dontshow{
write.ATetra(mygendata, file=tempfile())
}

}
\keyword{ file }
