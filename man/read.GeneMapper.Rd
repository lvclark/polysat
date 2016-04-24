\name{read.GeneMapper}
\alias{read.GeneMapper}
\title{
Read GeneMapper Genotypes Tables
}
\description{
  Given a vector of filepaths to tab-delimited text files containing
  genotype data in the ABI GeneMapper Genotypes Table format,
  \code{read.GeneMapper} produces a \code{genambig} object containing
  the genotype data.
}
\usage{
read.GeneMapper(infiles, forceInteger=TRUE)
}
\arguments{
  \item{infiles}{
A character vector of paths to the files to be read.
}
\item{forceInteger}{
  Boolean.  If \code{TRUE}, alleles will be coerced to integers.  This is
  particularly useful for stripping any white space from allels and
  preventing alleles from being imported as character strings.  If
  \code{FALSE}, alleles will be imported as numeric or character values,
  depending on the content of the input file(s).
  }
}
\details{
\code{read.GeneMapper} can read the genotypes tables that are exported
  by the Applied Biosystems GeneMapper software.  The only alterations
  to the files that the user may have to make are 1) delete
  any rows with missing data or fill in \code{-9} in the first allele slot for
  that row, 2) make sure that all allele names are numeric
  representations of fragment length (no question marks or dashes), and
  3) put sample names into the Sample Name column, if the names that you
  wish to use in analysis are not already there.  Each file should have
  the standard header row produced by the software.  If any sample has
  more than one genotype listed for a given locus, only the last
  genotype listed will be used.

  The file format is simple enough that the user can easily create files
  manually if GeneMapper is not the software used in allele calling.
  The files are tab-delimited text files.  There should be a header row
  with column names.  The column labeled \dQuote{Sample Name} should contain
  the names of the samples, and the column labeled \dQuote{Marker} should
  contain the names of the loci.  You can have as many or as few columns as
  needed to contain the alleles, and each of these columns should be
  labeled \dQuote{Allele X} where X is a number unique to each column.  Row
  labels and any other columns are ignored.  For any given sample, each
  allele is listed only once and is given as an integer that is the
  length of the fragment in nucleotides.  Alleles are separated by
  tabs.  If you have more allele columns than alleles for any given
  sample, leave the extra cells blank so that \code{read.table} will
  read them as \code{NA}.  Example data files in this format are
  included in the package.

  \code{read.GeneMapper} will read all of your data at once.  It takes
  as its first argument a character vector containing paths to all of
  the files to be read.  How the data are distributed over these files
  does not matter.  The function finds all unique sample names and all
  unique markers across all the files, and automatically puts a missing
  data symbol into the list if a particular sample and locus combination
  is not found.  Rows in which all allele cells are blank should NOT be
  included in the input files; either delete these rows or put the
  missing data symbol into the first allele cell.

  Sample and locus names must be consistent within and across the
  files.  The object that is produced is indexed by these names.

  If \code{forceInteger=FALSE}, alleles can be non-numeric values.  Some
  functionality of \pkg{polysat} will be lost in this case, but it could
  allow for the import of SNP data, for example.
}
\value{
A \code{genambig} object containing genotypes from the files, stored as
vectors of unique alleles in its \code{Genotypes} slot.  Other slots are
left at the default values.
}
\references{
%\url{http://www.appliedbiosystems.com/genemapper}
% Note that this URL redirects to a different page depending on the user's
% location.  Causes libcurl error code 52 in R CMD check.

GeneMapper website: \url{https://www.thermofisher.com/order/catalog/product/4475073}
}
\author{
Lindsay V. Clark
}
\note{
  A \sQuote{subscript out of bounds} error may mean that a sample name
  or marker was left blank in one of the input files.  A \sQuote{NAs
    introduced by coercion} warning when \code{forceInteger=TRUE} means that a
  non-numeric, non-whitespace character was included in one of the
  allele fields of the file(s), in which case the file(s) should be
  carefully checked and re-imported.
}

\seealso{
\code{\linkS4class{genambig}}, \code{\link{read.Structure}},
\code{\link{read.GenoDive}}, \code{\link{read.SPAGeDi}},
\code{\link{read.Tetrasat}}, \code{\link{read.ATetra}},
\code{\link{write.GeneMapper}}, \code{\link{read.POPDIST}},
\code{\link{read.STRand}}
}
\examples{
# create a table of data
gentable <- data.frame(Sample.Name=rep(c("ind1","ind2","ind3"),2),
                       Marker=rep(c("loc1","loc2"), each=3),
                       Allele.1=c(202,200,204,133,133,130),
                       Allele.2=c(206,202,208,136,142,136),
                       Allele.3=c(NA,208,212,145,148,NA),
                       Allele.4=c(NA,216,NA,151,157,NA)
                       )
# create a file (inspect this file in a text editor or spreadsheet
# software to see the required format)
write.table(gentable, file="readGMtest.txt", quote=FALSE, sep="\t",
            na="", row.names=FALSE, col.names=TRUE)

# read the file
mygenotypes <- read.GeneMapper("readGMtest.txt")

# inspect the results
viewGenotypes(mygenotypes)
}
\keyword{ file }

