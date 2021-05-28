\name{write.GeneMapper}
\alias{write.GeneMapper}
\title{
Write Genotypes to a Table Similarly to ABI GeneMapper
}
\description{
Given a \code{genambig} object,
\code{write.GeneMapper} writes a text file of a table containing columns
for sample name, locus, and alleles.
}
\usage{
write.GeneMapper(object, file = "", samples = Samples(object),
                 loci = Loci(object))
}
\arguments{
  \item{object}{
A \code{genambig} object containing genotype data to write to the file.
The \code{Ploidies} slot is used for determining how many allele
columns to make.
}
  \item{file}{
Character string.  The path to which to write the file.
}
  \item{samples}{
Character vector.  Samples to write to the file.  This should be a
subset of \code{Samples(object)}.
}
  \item{loci}{
Character vector.  Loci to write to the file.  This should be a subset
of \code{Loci(object)}.
}
}
\details{
  Although I do not know of any population genetic software other than
  \pkg{polysat} that will
  read this data format directly, the ABI GeneMapper Genotypes Table
  format is a convenient way for the user to store microsatellite
  genotype data and view it in a text editor or spreadsheet software.
  Each row contains the sample name, locus name, and alleles separated
  by tabs.

  The number of allele columns needed is detected by the maximum value
  of \code{Ploidies(object,samples,loci)}.  The function will add additional
  columns if it encounters genotypes with more than this number of
  alleles.

  \code{write.GeneMapper} handles missing data in a very simple way, in
  that it writes the missing data symbol directly to the table as though
  it were an allele.  If you want missing data to be represented
  differently in the table, you can open it in spreadsheet software and
  use find/replace or conditional formatting to locate missing data.

  The file that is produced can be read back into R directly by
  \code{\link{read.GeneMapper}}, and therefore may be a convenient way to backup
  genotype data for future analysis and manipulation in \pkg{polysat}.
  (\code{\link{save}} can also be used to backup an R object more directly,
  including population and other information.)  This can also enable the
  user to edit genotype data in spreadsheet software, if the
  \code{\link{editGenotypes}} function is not sufficient.

}
\value{
A file is written but no value is returned.
}
\references{
%  \url{http://www.appliedbiosystems.com/genemapper}
% Note that this URL redirects to a different page depending on the user's
% location.  Causes libcurl error code 52 in R CMD check.

GeneMapper website: \url{https://www.thermofisher.com/order/catalog/product/4475073}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.GeneMapper}}, \code{\link{write.Structure}},
\code{\link{write.GenoDive}}, \code{\link{write.Tetrasat}},
\code{\link{write.ATetra}}, \code{\link{write.POPDIST}},
\code{\link{write.SPAGeDi}}, \code{\link{editGenotypes}}
}
\examples{
# create a genotype object (usually done by reading a file)
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
Ploidies(mygendata) <- 6

\dontrun{
# write a GeneMapper file
write.GeneMapper(mygendata, "exampleGMoutput.txt")

# view the file with read.table
read.table("exampleGMoutput.txt", sep="\t", header=TRUE)
}
\dontshow{
write.GeneMapper(mygendata, tempfile())
}
}
\keyword{ file }
