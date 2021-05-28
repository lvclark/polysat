\name{write.GenoDive}
\alias{write.GenoDive}
\title{
Write a File in GenoDive Format
}
\description{
  \code{write.GenoDive} uses data from a \code{genambig} object to
  create a file formatted for the software GenoDive.
}
\usage{
write.GenoDive(object, digits = 2, file = "",
               samples = Samples(object), loci = Loci(object))
}
\arguments{
  \item{object}{
A \code{genambig} object containing genotypes, ploidies, population
identities, microsatellite repeat lengths, and description for the
dataset of interest.
}
  \item{digits}{
    An integer indicating how many digits to use to represent each allele
    (usually 2 or 3).
}
  \item{file}{
A character string of the file path to which to write.
}
  \item{samples}{
    A character vector of samples to include in the file.  A subset of
    \code{Samples(object)}.
}
  \item{loci}{
    A character vector of loci to include in the file.  A subset of
    \code{Loci(object)}.
}

}
\details{
The number of individuals, number of populations, number of loci, and
maximum ploidy of the sample are calculated automatically and entered in
the second line of the file.  If the maximum ploidy needs to be reduced
by random removal of alleles, it is possible to do this in the software
GenoDive after importing the data.  The \code{genambig} object should
not have individuals with more alleles than the highest ploidy level
listed in its \code{Ploidies} slot.

Several steps happen in order to convert alleles to the right format.
First, all instances of the missing data symbol are replaced with zero.
Alleles are then divided by the numbers provided in
\code{Usatnts(object)} (and rounded down if necessary) in order to
convert them from nucleotides to repeat numbers.  If the alleles are
still too big to be represented by the number of digits specified,
\code{write.GenoDive} repeatedly subtracts a number (\code{10^(digits-1)}; 10
if \code{digits=2}) from all alleles at a locus until the alleles are small
enough.  Alleles are then converted to characters, and a leading zero is
added to an allele if it does not have enough digits.  These alleles are
concatenated at each locus so that each sample*locus genotype is an
uninterrupted string of numbers.
}
\value{
A file is written but no value is returned.
}
\references{
  Meirmans, P. G. and Van Tienderen P. H. (2004) GENOTYPE and GENODIVE:
  two programs for the analysis of genetic diversity of asexual
  organisms. \emph{Molecular Ecology Notes} \bold{4}, 792-794.

  \url{http://www.bentleydrummer.nl/software/software/GenoDive.html}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.GenoDive}}, \code{\link{write.Structure}},
\code{\link{write.ATetra}}, \code{\link{write.Tetrasat}},
\code{\link{write.GeneMapper}}, \code{\link{write.POPDIST}},
\code{\link{write.SPAGeDi}}
}
\examples{
# set up the genotype object (usually done by reading a file)
mysamples <- c("Mal", "Inara", "Kaylee", "Simon", "River", "Zoe",
               "Wash", "Jayne", "Book")
myloci <- c("loc1", "loc2")
mygendata <- new("genambig", samples=mysamples, loci=myloci)
mygendata <- reformatPloidies(mygendata, output="sample")
Genotypes(mygendata, loci="loc1") <- list(c(304,306), c(302,310),
                                          c(306), c(312,314),
    c(312,314), c(308,310), c(312), c(302,308,310), c(-9))
Genotypes(mygendata, loci="loc2") <- list(c(118,133), c(121,130),
                                          c(122,139), c(124,133),
    c(118,124), c(121,127), c(124,136), c(124,127,136), c(121,130))
Usatnts(mygendata) <- c(2,3)
PopNames(mygendata) <- c("Core","Outer Rim")
PopInfo(mygendata) <- c(2,1,2,1,1,2,2,2,1)
Ploidies(mygendata) <- c(2,2,2,2,2,2,2,3,2)
Description(mygendata) <- "Serenity crew"

\dontrun{
# write files (use file="" to write to the console instead)
write.GenoDive(mygendata, digits=2, file="testGenoDive2.txt")
write.GenoDive(mygendata, digits=3, file="testGenoDive3.txt")
}
\dontshow{
write.GenoDive(mygendata, digits=2, file=tempfile())
write.GenoDive(mygendata, digits=3, file=tempfile())
}
}
\keyword{file}
