\name{read.POPDIST}
\alias{read.POPDIST}
\title{
Read Genotype Data in POPDIST Format
}
\description{
\code{read.POPDIST} reads one or more text files formatted for the
software POPDIST and produces a \code{"genambig"} object containing
genotypes, ploidies, and population identities from the file(s).
}
\usage{
read.POPDIST(infiles)
}
\arguments{
  \item{infiles}{
A character vector of file paths to be read.
}
}
\details{
  The format for the software POPDIST is a modified version of the
  popular Genepop format.  The first line is a comment line, followed by
  a list of locus names, each on a separate line or on one line
  separated by commas.  A line starting with
  the string \dQuote{Pop} (\dQuote{pop} and \dQuote{POP} are also recognized)
  indicates the beginning of data for one population.  Each individual
  is then represented on one line, with the population name and
  individual genotype separated by a tab followed by a comma.  Genotypes
  for different
  loci are separated by a tab or space.  Each allele must be coded by
  two digits.  Zeros (\dQuote{00}) indicate missing data, either for an
  entire locus or for a partially heterozygous genotype.  Partially
  heterozygous genotypes can also be represented by the arbitrary
  duplication of alleles.

If more than one file is read at once, locus names must be consistent
across all files.  Locus and population names should not start with \dQuote{Pop},
\dQuote{pop}, or \dQuote{POP}, as \code{read.POPDIST} searches for these
character strings in order to identify the lines that delimit populations.
}
\value{
A \code{"genambig"} object.  The \code{Description} slot of the object
is taken from the comment line of the first file.  Locus names are taken
from the files, and samples are given numbers instead of names.  Each
genotype consists of all unique non-zero integers for a given sample and
locus.  The \code{Ploidies} slot is filled in based on how many alleles
are present at each locus of each sample (the number of characters
for the genotype, divided by two).  \code{\link{reformatPloidies}} is
used internally by the function to collapse the ploidies to the simplest
format.  Population names are taken from the
individual genotype lines, and population identities are recorded based
on how the individuals are delimited by \dQuote{Pop} lines.
}
\references{
\url{http://gbi.agrsci.dk/~bernt/popgen/}

Tomiuk, J., Guldbrandtsen, B. and Loeschcke, B.  (2009)  Genetic
similarity of polyploids: a new version of the computer program POPDIST
(version 1.2.0) considers intraspecific genetic differentiation.
\emph{Molecular Ecology Resources} \bold{9}, 1364-1368.

Guldbrandtsen, B., Tomiuk, J. and Loeschcke, B.  (2000) POPDIST version
1.1.1: A program to calculate population genetic distance and identity
measures.  \emph{Journal of Heredity} \bold{91}, 178-179.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{write.POPDIST}}, \code{\link{read.Tetrasat}},
\code{\link{read.ATetra}}, \code{\link{read.Structure}},
\code{\link{read.SPAGeDi}}, \code{\link{read.GeneMapper}},
\code{\link{read.GenoDive}}, \code{\link{read.STRand}}
}
\examples{
# Create a file to read (this is typically done in a text editor)
cat("An example for the read.POPDIST documentation.",
"abcR",
"abcQ",
"Pop",
"Piscataqua\t, 0204 0505",
"Piscataqua\t, 0404 0307",
"Piscataqua\t, 050200 030509",
"Pop",
"Salmon Falls\t, 1006\t0805",
"Salmon Falls\t, 0510\t0308",
"Pop",
"Great Works\t, 050807 030800",
"Great Works\t, 0000 0408",
"Great Works\t, 0707 0305",
file="testPOPDIST.txt", sep="\n")

# View the file in the R console (or open it in a text editor)
cat(readLines("testPOPDIST.txt"), sep="\n")

# Read the file into a "genambig" object
fishes <- read.POPDIST("testPOPDIST.txt")

# View the data in the object
summary(fishes)
PopNames(fishes)
PopInfo(fishes)
Ploidies(fishes)
viewGenotypes(fishes)
}
\keyword{ file }
