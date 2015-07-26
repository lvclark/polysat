\name{write.SPAGeDi}
\alias{write.SPAGeDi}
\title{
Write Genotypes in SPAGeDi Format
}
\description{
  \code{write.SPAGeDi} uses data contained in a \code{genambig} object
  to create a file that can be read by the software SPAGeDi.  The
user controls how the genotypes are formatted, and can provide a data
frame of spatial coordinates for each sample.
}
\usage{
write.SPAGeDi(object, samples = Samples(object),
              loci = Loci(object), allelesep = "/",
              digits = 2, file = "",
              spatcoord = data.frame(X = rep(1, length(samples)),
                                     Y = rep(1, length(samples)),
                                     row.names = samples))
}
\arguments{
  \item{object}{
A \code{genambig} object containing genotypes, ploidies, population
identities, and microsatellite repeat lengths for the dataset of
interest.
}
  \item{samples}{
Character vector.  Samples to write to the file.  Must be a subset of
\code{Samples(object)}.
}
  \item{loci}{
Character vector.  Loci to write to the file.  Must be a subset of
\code{Loci(object)}.
}
  \item{allelesep}{
The character that will be used to separate alleles within a genotype.
If each allele should instead be a fixed number of digits, with no
characters to delimit alleles, set \code{allelesep = ""}.
}
  \item{digits}{
Integer.  The number of digits used to represent each allele.
}
  \item{file}{
A character string indicating the path to which the file should be
written.
}
  \item{spatcoord}{
Data frame.  Spatial coordinates of each sample.  Column names are used
for column names in the file.  Row names indicate sample, or if absent
it is assumed that the rows are in the same order as \code{samples}.
}
}
\details{
The Categories column of the SPAGeDi file that is produced contains
information from the \code{PopNames} and \code{PopInfo} slots of
\code{object}; the population name for each sample is written to the
column.

The first line of the file contains the number of individuals, number of
categories, number of spatial coordinates, number of loci, number of
digits for coding alleles, and maximum ploidy, and is generated
automatically from the data provided.

The function does not write distance intervals to the file, but instead
writes \code{0} to the second line.

All alleles for a given locus are divided by the \code{Usatnts} value
for that locus, after all missing data symbols have been replaced with
zeros.  If necessary, a multiple of 10 is subtracted from all
alleles at a locus in order to get the alleles down to the right number
of digits.

If a genotype has fewer alleles than the \code{Ploidies} value for
that sample and locus, zeros are added up to the ploidy.  If the
genotype has more
alleles than the ploidy, a random subset of alleles is used and a
warning is printed.  If the genotype has only one allele (is fully
heterozygous), then that allele is replicated to the ploidy of the
individual.  Genotypes are then concatenated into strings,
delimited by \code{allelesep}.  If \code{allelesep=""}, leading zeros
are first added to alleles as necessary to make them the right number of
digits.
}
\value{
A file is written but no value is returned.
}
\references{
  \url{http://ebe.ulb.ac.be/ebe/SPAGeDi.html}

  Hardy, O. J. and Vekemans, X. (2002) SPAGeDi: a versatile computer
  program to analyse spatial genetic structure at the individual or
  population levels.  \emph{Molecular Ecology Notes} \bold{2}, 618--620.
}
\author{
Lindsay V. Clark
}

\seealso{
  \code{\link{read.SPAGeDi}}, \code{\link{write.freq.SPAGeDi}},
  \code{\link{write.GenoDive}},
\code{\link{write.Structure}}, \code{\link{write.GeneMapper}},
\code{\link{write.ATetra}}, \code{\link{write.Tetrasat}}, 
\code{\link{write.POPDIST}}
}
\examples{
# set up data to write (usually read from a file)
mygendata <- new("genambig", samples = c("ind1","ind2","ind3","ind4"),
                 loci = c("loc1", "loc2"))
mygendata <- reformatPloidies(mygendata, output="sample")
Genotypes(mygendata, samples="ind1") <- list(c(102,106,108),c(207,210))
Genotypes(mygendata, samples="ind2") <- list(c(104),c(204,210))
Genotypes(mygendata, samples="ind3") <- list(c(100,102,108),c(201,213))
Genotypes(mygendata, samples="ind4") <- list(c(102,112),c(-9))
Ploidies(mygendata) <- c(3,2,2,2)
Usatnts(mygendata) <- c(2,3)
PopNames(mygendata) <- c("A", "B")
PopInfo(mygendata) <- c(1,1,2,2)
myspatcoord <- data.frame(X=c(27,29,24,30), Y=c(44,41,45,46),
                          row.names=c("ind1","ind2","ind3","ind4"))

# write a file
write.SPAGeDi(mygendata, spatcoord = myspatcoord,
              file="SpagOutExample.txt")
}
\keyword{ file }

