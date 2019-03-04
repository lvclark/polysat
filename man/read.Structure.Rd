\name{read.Structure}
\alias{read.Structure}
\title{
Read Genotypes and Other Data from a Structure File
}
\description{
\code{read.Structure} creates a \code{genambig} object by reading a text
file formatted for the software Structure.  \code{Ploidies} and
\code{PopInfo} (if
available) are also written to the object, and data from additional
columns can optionally be extracted as well.
}
\usage{
read.Structure(infile, ploidy, missingin = -9, sep = "\t",
               markernames = TRUE, labels = TRUE, extrarows = 1,
               popinfocol = 1, extracols = 1, getexcols = FALSE,
               ploidyoutput="one")
}
\arguments{
  \item{infile}{
Character string.  The file path to be read.
}
\item{ploidy}{
    Integer.  The ploidy of the file, \emph{i.e.} how many rows there
    are for each individual.
}
  \item{missingin}{
The symbol used to represent missing data in the Structure file.
}
  \item{sep}{
The character used to delimit the fields of the Structure file (tab by default).
}
  \item{markernames}{
Boolean, indicating whether the file has a header containing marker names.
}
  \item{labels}{
Boolean, indicating whether the file has a column containing sample names.
}
  \item{extrarows}{
Integer.  The number of extra rows that the file has, not counting
marker names.  This could include rows for recessive alleles,
inter-marker distances, or phase information.
}
\item{popinfocol}{
  Integer.  The column number (after the labels column, if present)
  where the data to be used for \code{PopInfo} are stored.  Can be
  \code{NA} to indicate that \code{PopInfo} should not be extracted from
  the file.
  }
  \item{extracols}{
Integer.  The number of extra columns that the file has, not counting
sample names (labels) but counting the column to be used for
\code{PopInfo}.  This could include PopData, PopFlag, LocData,
Phenotype, or any other extra columns.
}

  \item{getexcols}{
Boolean, indicating whether the function should return the data from any
extra columns.
}
\item{ploidyoutput}{
  This argument determines what is assigned to the \code{Ploidies} slot
  of the \code{"genambig"} dataset that is output.  It should be a
  string, either \code{"one"}, \code{"samplemax"}, or
  \code{"matrix"}.  This indicates, respectively, that \code{ploidy}
  should be used as the ploidy of the entire dataset, that the maximum
  number of alleles for each sample should be used as the ploidy of that
  sample, or that ploidies should be stored as a matrix of allele counts
  for each sample*locus.
  }
}
\details{
The current version of \code{read.Structure} does not support the
ONEROWPERIND option in the file format.  Each locus must only have one
column.  If your data are in ONEROWPERIND format, it should be fairly
simple to manipulate it in a spreadsheet program so that it can be read
by \code{read.GeneMapper} instead.

\code{read.Structure} uses \code{read.table} to initially read the file into a
data frame, then extracts information from the data frame.  Because of
this, any header rows (particularly the one containing marker names)
should have leading tabs (or spaces if \code{sep=" "}) so that the marker
names align correctly with their corresponding genotypes.  You should be
able to open the file in a spreadsheet program and have everything align
correctly.

If the file does not contain sample names, set \code{labels=FALSE}.  The
samples will be numbered instead, and if you like you can use the
\code{Samples<-} function to edit the sample names of the genotype object after
import.  Likewise, if \code{markernames=FALSE}, the
loci will be numbered automatically by the column names that
\code{read.table} creates, but these can also be edited after the fact.

The \code{Ploidies} slot of the \code{"genambig"} object that is created
is initially indexed by both sample and locus, with ploidy being
written to the slot on a per-genotype basis.  After all genotypes have
been imported, \code{\link{reformatPloidies}} is used to convert
\code{Ploidies} to the simplest possible format before the object is returned.
}
\value{
If \code{getexcols=FALSE}, the function returns only a \code{genambig} object.

If \code{getexcols=TRUE}, the function returns a list with two elements.  The
first, named \code{ExtraCol}, is a data frame, where the row names are the
sample names and each column is one of the extra columns from the file
(but with each sample only once instead of being repeated \code{ploidy}
number of times).  The second element is named \code{Dataset} and is the
genotype object described above.
}
\references{
  \url{https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf}

  Hubisz, M. J., Falush, D., Stephens, M. and Pritchard, J. K. (2009)
  Inferring weak population structure with the assistance of sample
  group information.  \emph{Molecular Ecology Resources} \bold{9},
  1322--1332.

  Falush, D., Stephens, M. and Pritchard, J. K.  (2007)
  Inferences of population structure using multilocus genotype data:
  dominant markers and null alleles.  \emph{Molecular Ecology Notes}
  \bold{7}, 574--578.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{write.Structure}}, \code{\link{read.GeneMapper}},
\code{\link{read.Tetrasat}}, \code{\link{read.ATetra}},
\code{\link{read.GenoDive}},
\code{\link{read.SPAGeDi}}, \code{\link{read.POPDIST}},
\code{\link{read.STRand}}
}
\examples{
# create a file to read (normally done in a text editor or spreadsheet
# software)
cat("\t\tRhCBA15\tRhCBA23\tRhCBA28\tRhCBA14\tRUB126\tRUB262\tRhCBA6\tRUB26",
    "\t\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "WIN1B\t1\t197\t98\t152\t170\t136\t208\t151\t99",
    "WIN1B\t1\t208\t106\t174\t180\t166\t208\t164\t99",
    "WIN1B\t1\t211\t98\t182\t187\t184\t208\t174\t99",
    "WIN1B\t1\t212\t98\t193\t170\t203\t208\t151\t99",
    "WIN1B\t1\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "WIN1B\t1\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "WIN1B\t1\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "WIN1B\t1\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD1\t2\t208\t100\t138\t160\t127\t202\t151\t124",
    "MCD1\t2\t208\t102\t153\t168\t138\t207\t151\t134",
    "MCD1\t2\t208\t106\t157\t180\t162\t211\t151\t137",
    "MCD1\t2\t208\t110\t159\t187\t127\t215\t151\t124",
    "MCD1\t2\t208\t114\t168\t160\t127\t224\t151\t124",
    "MCD1\t2\t208\t124\t193\t160\t127\t228\t151\t124",
    "MCD1\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD1\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD2\t2\t208\t98\t138\t160\t136\t202\t150\t120",
    "MCD2\t2\t208\t102\t144\t174\t145\t214\t150\t132",
    "MCD2\t2\t208\t105\t148\t178\t136\t217\t150\t135",
    "MCD2\t2\t208\t114\t151\t184\t136\t227\t150\t120",
    "MCD2\t2\t208\t98\t155\t160\t136\t202\t150\t120",
    "MCD2\t2\t208\t98\t157\t160\t136\t202\t150\t120",
    "MCD2\t2\t208\t98\t163\t160\t136\t202\t150\t120",
    "MCD2\t2\t208\t98\t138\t160\t136\t202\t150\t120",
    "MCD3\t2\t197\t100\t172\t170\t159\t213\t174\t134",
    "MCD3\t2\t197\t106\t174\t178\t193\t213\t176\t132",
    "MCD3\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD3\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD3\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD3\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD3\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    "MCD3\t2\t-9\t-9\t-9\t-9\t-9\t-9\t-9\t-9",
    sep="\n",file="structtest.txt")

# view the file
cat(readLines("structtest.txt"), sep="\n")

# read the structure file into genotypes and populations
testdata <- read.Structure("structtest.txt", ploidy=8)

# examine the results
testdata
}
\keyword{ file }
