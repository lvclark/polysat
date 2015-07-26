\name{genbinary-class}
\Rdversion{1.1}
\docType{class}
\alias{genbinary-class}
\alias{Absent,genbinary-method}
\alias{Absent<-,genbinary-method}
\alias{Genotype,genbinary-method}
\alias{Genotypes,genbinary-method}
\alias{Genotypes<-,genbinary-method}
\alias{initialize,genbinary-method}
\alias{Missing<-,genbinary-method}
\alias{Present,genbinary-method}
\alias{Present<-,genbinary-method}
\alias{Samples<-,genbinary-method}
\alias{Loci<-,genbinary-method}
\alias{isMissing,genbinary-method}
\alias{summary,genbinary-method}
\alias{editGenotypes,genbinary-method}
\alias{viewGenotypes,genbinary-method}
\alias{deleteSamples,genbinary-method}
\alias{deleteLoci,genbinary-method}
\alias{[,genbinary-method}
\alias{estimatePloidy,genbinary-method}

\title{Class "genbinary"}
\description{
This is a subclass of \code{gendata} that allows genotypes to be stored
as a matrix indicating the presence and absence of alleles.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("genbinary",
  samples, loci, ...)}.  After objects are initialized with sample and
locus names, data can be added to slots using the replacement functions.
}
\section{Slots}{
  \describe{
    \item{\code{Genotypes}:}{Object of class \code{"matrix"}.  Row names
    of the matrix are sample names.  Each column name is a locus name
    and an allele separated by a period (\emph{e.g.}
    \code{"loc1.124"}); each column represents on allele.  The number of
    alleles per locus is not limited and can be expanded even after
    entering initial data.  Each element of the matrix must be equal to
    either \code{Present(object)}, \code{Absent(object)}, or
    \code{Missing(object)}.  These symbols indicate, respectively, that
    a sample has an allele, that a sample does not have an allele, or
    that data for the sample at that locus are missing.}
    \item{\code{Present}:}{Object of class \code{"ANY"}.  The integer
      \code{1} by default.  This symbol is used in the \code{Genotypes}
      slot to indicate the presence of an allele in a sample. }
    \item{\code{Absent}:}{Object of class \code{"ANY"}.  The integer
      \code{0} by default.  This symbol is used in the \code{Genotypes}
      slot to indicate the absence of an allele in a sample. }
    \item{\code{Description}:}{Object of class \code{"character"}.  A
      character string or vector describing the dataset, for the
      convenience of the user. }
    \item{\code{Missing}:}{Object of class \code{"ANY"}.  The integer
      \code{-9} by default.  This symbol is used in the \code{Genotypes}
    slot to indicate that data are missing for a given sample and locus.}
    \item{\code{Usatnts}:}{Object of class \code{"integer"}.  A vector,
      named by loci.  This indicates the repeat length of each locus.
      \code{2} indicates dinucleotide repeats, \code{3}
      indicates trinucleotide repeats, and so on.  If the alleles stored
      in the column names of the \code{Genotypes} slot for a given locus
      are already written
    in terms of repeat number, the \code{Usatnts} value for that locus
    should be \code{1}.  In other words, all alleles for a locus can be
    divided by the number in \code{Usatnts} to give alleles expressed in
  terms of relative repeat number. }
    \item{\code{Ploidies}:}{Object of class \code{"integer"}.  A vector,
    named by samples.  This indicates the ploidy of each sample.}
    \item{\code{PopInfo}:}{Object of class \code{"integer"}.  A vector,
      named by samples.  This indicates the population identity of each
      sample. }
    \item{\code{PopNames}:}{Object of class \code{"character"}.  Names
      of each population.  The position of the population name in the
      vector corresponds to the number used to represent that population
      in the \code{PopInfo} slot.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{gendata}"}, directly.
}
\section{Methods}{
  \describe{
    \item{Absent}{\code{signature(object = "genbinary")}: Returns the
      symbol used to indicate that a given allele is absent in a given
      sample. }
    \item{Absent<-}{\code{signature(object = "genbinary")}: Changes the
      symbol used to indicate that a given allele is absent in a given
      sample.  The matrix in the \code{Genotypes} slot is searched for
      the old symbol, which is replaced by the new.  The new symbol is
      then written to the \code{Absent} slot. }
    \item{Genotype}{\code{signature(object = "genbinary")}: Returns a
      matrix containing the genotype for a given sample and locus (by a
      call to \code{Genotypes}). }
    \item{Genotypes}{\code{signature(object = "genbinary")}: Returns the
    matrix stored in the \code{Genotypes} slot, or a subset as specified
  by the \code{samples} and \code{loci} arguments.}
    \item{Genotypes<-}{\code{signature(object = "genbinary")}: A method
      for adding or replacing genotype data in the object.  Note that
      allele columns cannot be removed from the matrix in the
      \code{Genotypes} slot using this method, although an entire column
      could be filled with zeros in order to effectively remove an
      allele from the dataset.  If the order of rows in \code{value}
      (the matrix containing values to be assigned to the
      \code{Genotypes} slot) is not identical to \code{Samples(object)},
    the \code{samples} argument should be used to indicate row order.
    Row names in \code{value} are ignored.  The \code{loci} argument can
  be left at the default, even if only a subset of loci are being
  assigned.  Column names of \code{value} are important, and should be
  the locus name and allele name separated by a period, as they are in
  the \code{Genotypes} slot.  After checking that the column name is
  valid, the method checks for whether the column name already exists or
not in the \code{Genotypes} slot.  If it does exist, data from that
column are replaced with data from \code{value}.  If not, a column is
added to the matrix in the \code{Genotypes} slot for the new allele.  If
the column is new and data are not being written for samples, the method
automatically fills in \code{Missing} or \code{Absent} symbols for
additional samples, depending on whether or not data for the locus
appear to be missing for the sample or not.}
    \item{initialize}{\code{signature(.Object = "genbinary")}: Sets up a
    \code{genbinary} object when \code{new("genbinary")} is called.  If
    \code{samples} or \code{loci} arguments are missing, these are
    filled in with dummy values (\code{"ind1", "ind2", "loc1",
      "loc2"}).  The matrix is then set up in the \code{Genotypes}
    slot.  Sample names are used for row names, and there are zero
    columns.  The \code{initialize} method for \code{gendata} is then
    called.}
    \item{Missing<-}{\code{signature(object = "genbinary")}: Replaces
      all elements in matrix in the \code{Genotypes} slot containing the
    old \code{Missing} symbol with the new \code{Missing} symbol.  The
    method for \code{gendata} is then called to replace the value in the
  \code{Missing} slot.}
    \item{Present}{\code{signature(object = "genbinary")}: Returns the
      symbol used to indicate that a given allele is present in a given
      sample. }
    \item{Present<-}{\code{signature(object = "genbinary")}: Changes the
    symbol used for indicating that a given allele is present in a given
  sample.  The symbol is first replaced in the \code{Genotypes} slot,
  and then in the \code{Present} slot.}
\item{Samples<-}{\code{signature(object = "genbinary")}: Changes sample
  names in the dataset.  Changes the row names in the \code{Genotypes}
  slot, then calls the method for \code{gendata} to change the names in
  the \code{PopInfo} and \code{Ploidies} slots.}
\item{Loci<-}{\code{signature(object = "genbinary")}: Changes locus
  names in the dataset.  Replaces the locus portion of the column names
  in the \code{Genotypes} slot, then calls the method for \code{gendata}
  to change the names in the \code{Usatnts} slot.}
\item{isMissing}{\code{signature(object = "genbinary")}: Returns Boolean
values, by sample and locus, indicating whether genotypes are missing.
If there are any missing data symbols within the genotype, it is
considered missing.}
\item{summary}{\code{signature(object = "genbinary")}: Prints
  description of dataset and number of missing genotypes, then calls the
  method for \code{gendata} to print additional information.}
\item{editGenotypes}{\code{signature(object = "genbinary")}: Opens the
  genotype matrix in the Data Editor for editing.  Useful for making
  minor changes, although allele columns cannot be added using this method.}
\item{viewGenotypes}{\code{signature(object = "genbinary")}: Prints the
  genotype matrix to the console, one locus at a time.}
\item{deleteSamples}{\code{signature(object = "genbinary")}: Removes the
specified samples from the genotypes matrix, then calls the method for
\code{gendata}.}
\item{deleteLoci}{\code{signature(object = "genbinary")}: Removes the
  specified loci from the genotypes matrix, then calls the method for
  \code{gendata}.}
\item{"["}{\code{signature(x = "genbinary", i = "ANY", j = "ANY")}:
  Subscripting method.  Returns a \code{genbinary} object with a subset
  of the samples and/or loci from \code{x}.  Usage:
  \code{genobject[samples,loci]}.}
\item{estimatePloidy}{\code{signature(object = "genbinary")}: Creates a
  data frame of mean and maximum number of alleles per sample, which is
  then opened in the Data Editor so that the user can manually specify
  the ploidy of each sample.  Ploidies are then written to the
  \code{Ploidies} slot of the object.}
\item{merge}{\code{signature(x = "genbinary", y = "genbinary")}: Merges
  two genotype objects together.  See
  \code{\link{merge,genbinary,genbinary-method}}.}
	 }
}

\author{
Lindsay V. Clark
}


\seealso{
\code{\linkS4class{gendata}}, \code{\link[=Samples]{Accessors}},
\code{\linkS4class{genambig}}
}
\examples{
# show the class definition
showClass("genbinary")

# create a genbinary object
mygen <- new("genbinary", samples = c("indA", "indB", "indC", "indD"),
             loci = c("loc1", "loc2"))
Description(mygen) <- "Example genbinary object for the documentation."
Usatnts(mygen) <- c(2,3)
PopNames(mygen) <- c("Maine", "Indiana")
PopInfo(mygen) <- c(1,1,2,2)
Genotypes(mygen) <- matrix(c(1,1,0,0, 1,0,0,1, 0,0,1,1,
                             1,-9,1,0, 0,-9,0,1, 1,-9,0,1, 0,-9,1,1),
   nrow=4, ncol=7, dimnames = list(NULL,
   c("loc1.140", "loc1.144", "loc1.150",
     "loc2.97", "loc2.100", "loc2.106", "loc2.109")))

# view all of the data in the object
mygen
}
\keyword{classes}
