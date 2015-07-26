\name{write.POPDIST}
\alias{write.POPDIST}
\title{
Write Genotypes to a POPDIST File
}
\description{
\code{write.POPDIST} uses data from a \code{"genambig"} object to write
a file formatted for the software POPDIST.
}
\usage{
write.POPDIST(object, samples = Samples(object),
              loci = Loci(object), file = "")
}
\arguments{
  \item{object}{
A \code{"genambig"} object.  \code{Ploidies} and \code{PopInfo} are
required, and if provided \code{Usatnts} may be used to convert alleles
to repeat number in order to represent each allele with two digits.
Locus names and \code{PopNames} are used in the file, but sample
names are not.
}
  \item{samples}{
An optional character vector of samples to use.  Must be a subset of
\code{Samples(object)}.
}
  \item{loci}{
An optional character vector of loci to use.  Must be a subset of
\code{Loci(object)}.
}
  \item{file}{
Character string.  File path to which to write.
}
}
\details{
POPDIST is a program that calculates inter-population distance measures,
some of which are available for polyploid samples with allele copy
number ambiguity.  Each population must be of uniform ploidy, but
different populations may have different ploidies.

Two types of warning messages may be printed by \code{write.POPDIST}.
The first indicates that a population contains individuals of more than
one ploidy.  In this case a file is still written, but POPDIST may not
be able to read it.  Separate populations with different ploidies are
okay.  The second type of warning indicates that an individual has more
alleles than its ploidy level.  If this occurs, alleles are randomly
removed from the genotype that is written to the file.

If necessary, \code{write.POPDIST} converts alleles into a two-digit
format, similarly to \code{write.Tetrasat}.  If the value of any allele
for a given locus is greater than 99, the function first checks to see
if the locus has a \code{Usatnts} value greater than 1, and if so
divides all alleles by this value and rounds down.  If the locus still
has alleles with more than two digits, a multiple of 10 is subtracted
from all alleles.  A zero is placed in front of any allele with one digit.
}
\value{
A file is written but no value is returned.
}
\references{
  \url{http://gbi.agrsci.dk/~bernt/popgen/}

Tomiuk, J., Guldbrandtsen, B. and Loeschcke, B.  (2009)  Genetic
similarity of polyploids: a new version of the computer program POPDIST
(version 1.2.0) considers intraspecific genetic differentiation.
\emph{Molecular Ecology Resources} \bold{9}, 1364-1368.

Guldbrandtsen, B., Tomiuk, J. and Loeschcke, B.  (2000) POPDIST version
1.1.1: A program to calculate population genetic distance and identity
measures.  \emph{Journal of Heredity} \bold{91}, 178--179.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.POPDIST}}, \code{\link{write.Tetrasat}},
\code{\link{write.ATetra}}, \code{\link{write.SPAGeDi}},
\code{\link{write.GenoDive}}, \code{\link{write.Structure}},
\code{\link{write.GeneMapper}}
}
\examples{
# create a "genambig" object containing the dataset
mygen <- new("genambig", samples=c("a", "b", "c", "d"),
             loci=c("loc1", "loc27"))
mygen <- reformatPloidies(mygen, output="sample")
Description(mygen) <- "Some example data for POPDIST"
PopInfo(mygen) <- c(1,1,2,2)
PopNames(mygen) <- c("Old Orchard Beach", "York Beach")
Ploidies(mygen) <- c(2,2,4,4)
Usatnts(mygen) <- c(2,2)
Genotypes(mygen, loci="loc1") <- list(c(128, 134), c(130),
                                      Missing(mygen), c(126, 128, 132))
Genotypes(mygen, loci="loc27") <- list(c(209,211), c(207,217),
                                       c(207,209,215,221), c(211,223))

# write the file
write.POPDIST(mygen, file="forPOPDIST.txt")

# view the file
cat(readLines("forPOPDIST.txt"), sep="\n")

}
\keyword{ file }
