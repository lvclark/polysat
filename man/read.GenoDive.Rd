\name{read.GenoDive}
\alias{read.GenoDive}
\title{
Import Genotype Data from GenoDive File
}
\description{
\code{read.GenoDive} takes a text file in the format for the software
GenoDive and produces a \code{genambig} object.
}
\usage{
read.GenoDive(infile)
}
\arguments{
  \item{infile}{
A character string.  The path to the file to be read.
}
}
\details{
GenoDive is a Mac-only program for population genetic analysis that
allows for polyploid data.  \code{read.GenoDive} imports data from text
files formatted for this program.

The first line of the file is a comment line, which is written to the
\code{Description} slot of the \code{genambig} object.
On the second line, separated by tabs, are the number of individuals,
number of populations, number of loci, maximum ploidy
(ignored), and number of digits used to code alleles.

The following lines contain the names of populations, which are written
to the \code{PopNames} slot of the \code{genambig} object.  After that is a
header line for the genotype data.  This line contains, separated by
tabs, column headers for populations, clones (optional), and
individuals, followed by the name of each locus.  The locus names for the
genotype object are derived from this line.

Each individual is on one line following the genotype header line.
Separated by tabs are the population number, the clone number
(optional), the individual name (used as the sample name in the output)
and the genotypes at each locus.  Alleles at one locus are concatenated
together in one string without any characters to separate them.  Each
allele must have the same number of digits, although leading zeros can
be omitted.

If the only alleles listed for a particular individual and locus are
zeros, this is interpreted by \code{read.GenoDive} as missing data, and
\code{Missing(object)} (the default, \code{-9}) is written in that
genotype slot in the \code{genambig} object.
GenoDive allows for a genotype to be partially missing but \pkg{polysat} does
not; therefore, if an allele is coded as zero but other alleles are
recorded for that sample and locus, the output genotype will just
contain the alleles that are present, with the zeros thrown out.
}
\value{
A \code{genambig} object containing the data from the file.
}
\references{
  Meirmans, P. G. and Van Tienderen, P. H. (2004) GENOTYPE and GENODIVE:
  two programs for the analysis of genetic diversity of asexual
  organisms. \emph{Molecular Ecology Notes} \bold{4}, 792-794.

  \url{http://www.bentleydrummer.nl/software/software/GenoDive.html}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.GeneMapper}}, \code{\link{write.GenoDive}},
\code{\link{read.Tetrasat}}, \code{\link{read.ATetra}},
\code{\link{read.Structure}},
\code{\link{read.SPAGeDi}}, \code{\link{read.POPDIST}},
\code{\link{read.STRand}}
}
\examples{
# create data file (normally done in a text editor or spreadsheet software)
myfile <- tempfile()
cat(c("example comment line", "5\t2\t2\t3\t2", "pop1", "pop2",
      "pop\tind\tloc1\tloc2", "1\tJohn\t102\t1214",
      "1\tPaul\t202\t0", "2\tGeorge\t101\t121213",
      "2\tRingo\t10304\t131414","1\tYoko\t10303\t120014"),
    file = myfile, sep = "\n")

# import file data
exampledata <- read.GenoDive(myfile)

# view data
summary(exampledata)
viewGenotypes(exampledata)
exampledata
}
\keyword{file}
