\name{write.Tetrasat}
\alias{write.Tetrasat}
\title{
Write Genotype Data in Tetrasat Format
}
\description{
  Given a \code{genambig} object,
  \code{write.Tetrasat} creates a file that can be read by the software
  Tetrasat and Tetra.
}
\usage{
write.Tetrasat(object, samples = Samples(object),
               loci = Loci(object), file = "")
}
\arguments{
  \item{object}{
A \code{genambig} object containing the dataset of interest.  Genotypes,
population identities, microsatellite repeat lengths, and the dataset
description of \code{object} are used by the function.
}
  \item{samples}{
A character vector of samples to write to the file.  Should be a subset
of \code{Samples(object)}.
}
  \item{loci}{
A character vector of loci to write to the file.  Should be a subset of
\code{Loci(object)}.
}
  \item{file}{
A character string indicating the file to which to write.
}
}
\details{
Tetrasat files are space-delimited text files in which all alleles at a
locus are concatenated into a string eight characters long.  Population
names or numbers are not used in the file, but samples are ordered by
population, with the line \dQuote{Pop} delimiting populations.

\code{write.Tetrasat} divides each allele by the length of the repeat and
rounds down in order to convert alleles to repeat numbers.  If
necessary, it subtracts a multiple of 10 from all alleles at a locus to
make all allele values less than 100, or puts a zero in front of the
number if it only has one digit.  If the individual is fully homozygous
at a locus, the single allele is repeated four times.  If any genotype
has more than four alleles, \code{write.Tetrasat} picks a random sample of
four alleles without replacement, and prints a warning.  Missing data
are represented by blank spaces.

Sample names should be a maximum of 20 characters long in order for the
file to be read correctly by Tetrasat or Tetra.
}
\value{
A file is written but no value is returned.
}
\references{
  Markwith, S. H., Stewart, D. J. and Dyer, J. L.  (2006)
  TETRASAT: a program for the population analysis of
  allotetraploid microsatellite data. \emph{Molecular Ecology Notes}
  \bold{6}, 586-589.

  \url{http://ecology.bnu.edu.cn/zhangdy/TETRA/TETRA.htm}

  Liao, W. J., Zhu, B. R., Zeng, Y. F. and Zhang, D. Y.  (2008)
  TETRA: an improved program for population genetic analysis of
  allotetraploid microsatellite data.  \emph{Molecular Ecology Resources}
  \bold{8}, 1260--1262.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.Tetrasat}},
\code{\link{write.GeneMapper}}, \code{\link{write.ATetra}},
\code{\link{write.POPDIST}}
}
\examples{
# set up sample data (usually done by reading files)
mysamples <- c("ind1", "ind2", "ind3", "ind4")
myloci <- c("loc1", "loc2")
mygendata <- new("genambig", samples = mysamples, loci = myloci)
mygendata <- reformatPloidies(mygendata, output="one")
Usatnts(mygendata) <- c(2, 3)
Genotypes(mygendata, loci="loc1") <- list(c(202,204), c(204),
                                          c(200,206,208,212),
                                          c(198,204,208))
Genotypes(mygendata, loci="loc2") <- list(c(78,81,84),
                                          c(75,90,93,96,99),
                                          c(87), c(-9))
PopInfo(mygendata) <- c(1,2,1,2)
Description(mygendata) <- "An example for write.Tetrasat."
Ploidies(mygendata) <- 4

# write a Tetrasat file
write.Tetrasat(mygendata, file="tetrasattest.txt")

# view the file
cat(readLines("tetrasattest.txt"),sep="\n")
}
\keyword{file}
