\name{read.SPAGeDi}
\alias{read.SPAGeDi}
\title{
Read Genotypes in SPAGeDi Format
}
\description{
\code{read.SPAGeDi} can read a text file formatted for the SPAGeDi
software and return a \code{genambig} object,
as well as optionally returning a data frame of spatial coordinates.
The \code{genambig} object includes genotypes, ploidies, and population
identities (from the category column, if present) from the file.
}
\usage{
read.SPAGeDi(infile, allelesep = "/", returnspatcoord = FALSE)
}
\arguments{
  \item{infile}{
A character string indicating the path of the file to read.
}
  \item{allelesep}{
The character that is used to delimit alleles within genotypes, or
\code{""} if alleles have a fixed number of digits and are not delimited
by any character.  Other examples shown in section 3.2.1 of the SPAGeDi
1.3 manual include \code{"/"}, \code{" "}, \code{", "}, \code{"."}, and
\code{"--"}.
}
  \item{returnspatcoord}{
Boolean.  Indicates whether a data frame should be returned containing
the spatial coordinates columns.
}
}
\details{
SPAGeDi offers a lot of flexibility in how data files are formatted.
\code{read.SPAGeDi} accomodates most of that flexibility.  The primary
exception is that alleles must be delimited in the same way across all
genotypes, as specified by \code{allelesep}.  Comment lines beginning
with \code{//}, as well as blank lines, are ignored by
\code{read.SPAGeDi} just as they are by SPAGeDi.

\code{read.SPAGeDi} is not designed to read dominant data (see section
3.2.2 of the SPAGeDi 1.3 manual).  However, see
\code{genbinary.to.genambig} for a way to read this type
of data after some simple manipulation in a spreadsheet program.

The first line of a SPAGeDi file contains information that is used by
\code{read.SPAGeDi}.  The ploidy as specified in the 6th position of the
first line is ignored, and is instead calculated by counting alleles for
each individual (including zeros on the right, but not the left, side of
the genotype).  The number of digits specified in the 5th position of
the first line is only used if \code{allelesep=""}.  All other values
in the first line are important for the function.

If the only alleles found for a particular individual and locus are
zeros, the genotype is interpreted as missing.  Otherwise, zeros on the
left side of a genotype are ignored, and zeros on the right side of a
genotype are used in calculating the ploidy but are not included in the
genotype object that is returned.  If \code{allelesep=""},
\code{read.SPAGeDi} checks that the number of characters in the genotype
can be evenly divided by the number of digits per allele.  If not, zeros
are added to the left of the genotype string before splitting it into
alleles.

The \code{Ploidies} slot of the \code{"genambig"} object that is created
is initially indexed by both sample and locus, with ploidy being
written to the slot on a per-genotype basis.  After all genotypes have
been imported, \code{\link{reformatPloidies}} is used to convert
\code{Ploidies} to the simplest possible format before the object is returned.
}
\value{
  Under the default where \code{returnspatcoord=FALSE}, a \code{genambig}
  object is returned.  Alleles are formatted as integers.  The
  \code{Ploidies} slot is filled in according to the number of alleles
  per genotype, ignoring zeros on the left.  If the first line of the
  file indicates that there are more than zero categories, the category
  column is used to fill in the \code{PopNames} and \code{PopInfo} slots.

Otherwise, a list is returned:
\item{SpatCoord }{A data frame of spatial coordinates,
unchanged from the file.  The format of each column is determined under
the default \code{read.table} settings.  Row names are individual names
from the file.  Column names are the same as in the file.}
\item{Dataset }{A \code{genambig} object as described above.}
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
\code{\link{write.SPAGeDi}}, \code{\link{genbinary.to.genambig}},
\code{\link{read.table}}, \code{\link{read.GeneMapper}},
\code{\link{read.GenoDive}}, \code{\link{read.Structure}},
\code{\link{read.ATetra}}, \code{\link{read.Tetrasat}},
\code{\link{read.POPDIST}}, \code{\link{read.STRand}}
}
\examples{
# create a file to read (usually done with spreadsheet software or a
# text editor):
cat("// here's a comment line at the beginning of the file",
"5\t0\t-2\t2\t2\t4",
"4\t5\t10\t50\t100",
"Ind\tLat\tLong\tloc1\tloc2",
"ind1\t39.5\t-120.8\t00003133\t00004040",
"ind2\t39.5\t-120.8\t3537\t4246",
"ind3\t42.6\t-121.1\t5083332\t40414500",
"ind4\t38.2\t-120.3\t00000000\t41430000",
"ind5\t38.2\t-120.3\t00053137\t00414200",
"END",
sep="\n", file="SpagInputExample.txt")

# display the file
cat(readLines("SpagInputExample.txt"), sep="\n")

# read the file
mydata <- read.SPAGeDi("SpagInputExample.txt", allelesep = "",
returnspatcoord = TRUE)

# view the data
mydata
viewGenotypes(mydata[[2]])

}
\keyword{ file }

