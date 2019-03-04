\name{write.Structure}
\alias{write.Structure}
\title{Write Genotypes in Structure 2.3 Format}
\description{
  Given a dataset stored in a \code{genambig} object,
  \code{write.Structure} produces a text file of the genotypes in a format
  readable by Structure 2.2 and higher.  The user specifies the overall
  ploidy of the file, while the ploidy of each sample is extracted from
  the \code{genambig} object.  \code{PopInfo} and other data can
  optionally be written to the file as well.
}
\usage{
write.Structure(object, ploidy, file="",
                samples = Samples(object), loci = Loci(object),
                writepopinfo = TRUE, extracols = NULL,
                missingout = -9)
}
\arguments{
  \item{object}{A \code{genambig} object containing the data to write to
    the file.  There must be non-\code{NA} values of \code{Ploidies}
    (and \code{PopInfo} if \code{writepopinfo == TRUE}) for
    \code{samples}.}
  \item{ploidy}{PLOIDY for Structure, \emph{i.e.} how many rows per individual
    to write.}
  \item{file}{A character string specifying where the file should be written.}
  \item{samples}{An optional character vector listing the names of samples to be
    written to the file.}
  \item{loci}{An optional character vector listing the names of the loci to be
    written to the file.}
  \item{writepopinfo}{\code{TRUE} or \code{FALSE}, indicating whether to
  write values from the \code{PopInfo} slot of \code{object} to the file.}
  \item{extracols}{An array, with the first dimension names
    corresponding to \code{samples}, of PopData, PopFlag, LocData, Phenotype,
    or other values to be included in the extra columns in the file.}
  \item{missingout}{The number used to indicate missing data.}
}
\value{
  No value is returned, but instead a file is written at the path specified.
}
\details{
  Structure 2.2 and higher can process autopolyploid microsatellite data,
  although 2.3.3 or higher is recommended for its improvements on
  polyploid handling.  The input format of Structure requires that
  each locus take up one column and that each individual take up as
  many rows as the parameter PLOIDY.  Because of the multiple rows per
  sample, each sample name must be duplicated, as well as any
  population, location, or phenotype data.  Partially heterozygous
  genotypes also must have one arbitrary allele duplicated up to the
  ploidy of the sample, and samples that have a lower ploidy than that
  used in the file (for mixed polyploid data sets) must have a missing
  data symbol inserted to fill in the extra rows.  Additionally, if
  some samples have more alleles than PLOIDY (if you are using a lower
  PLOIDY to save processing time, or if there are extra alleles from
  scoring errors), some alleles must be randomly removed from the data.
  \code{write.Structure} performs this duplication, insertion, and random
  deletion of data.

The sample names from \code{samples} will be used as row
  names in the Structure file.  Each sample name should only be in the
  vector \code{samples} once, because \code{write.Structure} will duplicate
  the sample names a number of times as dictated by \code{ploidy}.

  In writing genotypes to the file, \code{write.Structure} compares the number
  of alleles in the genotype, the ploidy of the sample*locus as stored in
  \code{Ploidies}, and the ploidy of the file as stored in
  \code{ploidy}, and does one of six things (for a given sample x and
  locus loc):

  1) If \code{Ploidies(object,x,loc)} is greater than or equal to
  \code{ploidy}, and
  \code{length(Genotype(object, x, loc))} is equal to \code{ploidy}, the
  genotype data are used as is.

  2) If \code{Ploidies(object,x,loc)} is greater than or equal to
  \code{ploidy}, and
  \code{length(Genotype(object, x, loc))} is less than \code{ploidy},
  the first allele is
  duplicated as many times as necessary for there to be as many alleles
  as \code{ploidy}.

  3) If \code{Ploidies(object,x,loc)} is greater than or equal to
  \code{ploidy}, and \code{length(Genotype(object, x, loc))} is greater
  than \code{ploidy}, a random sample of
  the alleles, without replacement, is used as the genotype.

  4) If \code{Ploidies(object,x,loc)} is less than \code{ploidy}, and
  \code{length(Genotype(object, x, loc))} is equal to
  \code{Ploidies(object,x,loc)}, the genotype
  data are used as is and missing data symbols are inserted in the extra
  rows.

  5) If \code{Ploidies(object,x,loc)} is less than \code{ploidy}, and
  \code{length(Genotype(object, x, loc))} is less than
  \code{Ploidies(object,x,loc)}, the first
  allele is duplicated as many times as necessary for there to be as
  many alleles as \code{Ploidies(object,x,loc)}, and missing data symbols
  are inserted in the extra rows.

  6) If \code{Ploidies(object,x,loc)} is less than \code{ploidy}, and
  \code{length(Genotype(object, x, loc))} is greater than
  \code{Ploidies(object,x,loc)}, a random
  sample, without replacement, of \code{Ploidies(object)[x]} alleles is
  used, and
  missing data symbols are inserted in the extra rows. (Alleles are
  removed even though there is room for them in the file.)

  Two of the header rows that are optional for Structure are written by
  \code{write.Structure}.  These are \sQuote{Marker Names}, containing the
  names of loci supplied in \code{gendata}, and \sQuote{Recessive
    Alleles}, which contains the missing data symbol once for each locus.
  This indicates to the program that all alleles are codominant with
  copy number ambiguity.
}
\references{
  \url{https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf}

  Hubisz, M. J., Falush, D., Stephens, M. and Pritchard, J. K. (2009)
  Inferring weak population structure with the assistance of sample
  group information.  \emph{Molecular Ecology Resources} \bold{9}, 1322-1332.

  Falush, D., Stephens, M. and Pritchard, J. K.  (2007)
  Inferences of population structure using multilocus genotype data:
  dominant markers and null alleles.  \emph{Molecular Ecology Notes}
  \bold{7}, 574-578.
}
\note{If \code{extracols} is a character array, make sure none of the
  elements contain white space.}
\seealso{
  \code{\link{read.Structure}}, \code{\link{write.GeneMapper}},
  \code{\link{write.GenoDive}}, \code{\link{write.SPAGeDi}},
  \code{\link{write.ATetra}}, \code{\link{write.Tetrasat}},
\code{\link{write.POPDIST}}

}
\examples{
# input genotype data (this is usually done by reading a file)
mygendata <- new("genambig", samples = c("ind1","ind2","ind3",
                                         "ind4","ind5","ind6"),
                 loci = c("locus1","locus2"))
Genotypes(mygendata) <- array(list(c(100,102,106,108,114,118),c(102,110),
                      c(98,100,104,108,110,112,116),c(102,106,112,118),
                      c(104,108,110),c(-9),
                      c(204),c(206,208,210,212,220,224,226),
                      c(202,206,208,212,214,218),c(200,204,206,208,212),
                      c(-9),c(202,206)),
                 dim=c(6,2))
Ploidies(mygendata) <- c(6,6,6,4,4,4)
# Note that some of the above genotypes have more or fewer alleles than
# the ploidy of the sample.

# create a vector of sample names to be used.  Note that this excludes
#  ind6.
mysamples <- c("ind1","ind2","ind3","ind4","ind5")

# Create an array containing data for additional columns to be written
# to the file.  You might also prefer to just read this and the ploidies
# in from a file.
myexcols <- array(data=c(1,2,1,2,1,1,1,0,0,0),dim=c(5,2),
                  dimnames=list(mysamples, c("PopData","PopFlag")))

# Write the Structure file, with six rows per individual.
# Since outfile="", the data will be written to the console instead of a file.
write.Structure(mygendata, 6, "", samples = mysamples, writepopinfo = FALSE,
                extracols = myexcols)
}
\author{Lindsay V. Clark}
\keyword{file}