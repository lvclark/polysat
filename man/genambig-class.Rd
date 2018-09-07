\name{genambig-class}
\Rdversion{1.1}
\docType{class}
\alias{genambig-class}
\alias{deleteLoci,genambig-method}
\alias{deleteSamples,genambig-method}
\alias{editGenotypes,genambig-method}
\alias{estimatePloidy,genambig-method}
\alias{Genotype,genambig-method}
\alias{Genotype<-,genambig-method}
\alias{Genotypes,genambig-method}
\alias{Genotypes<-,genambig-method}
\alias{initialize,genambig-method}
\alias{isMissing,genambig-method}
\alias{Loci<-,genambig-method}
\alias{Missing<-,genambig-method}
\alias{Samples<-,genambig-method}
\alias{summary,genambig-method}
\alias{show,genambig-method}
\alias{viewGenotypes,genambig-method}
\alias{[,genambig-method}

\title{Class "genambig"}
\description{
Objects of this class store microsatellite datasets in which allele copy
number is ambiguous.  Genotypes are stored as a two-dimensional list of
vectors, each vector containing all unique alleles for a given sample at
a given locus.  \code{genambig} is a subclass of \code{gendata}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("genambig",
  samples, loci, ...)}.  This automatically sets up a two-dimensional
list in the \code{Genotypes} slot, with \code{dimnames=list(samples,
  loci)}.  This array-list is initially populated with the missing data
symbol.  All other slots are given initial values according to the
\code{initialize} method for \code{gendata}.  Data can then be inserted
into the slots using the replacement functions (see
\code{\link[=Samples]{Accessors}}).

}
\section{Slots}{
  \describe{
    \item{\code{Genotypes}:}{Object of class \code{"array"}.  The first
      dimension of the array represents and is named by samples, while
      the second dimension represents and is named by loci.  Each
      element of the array can contain a vector.  Each vector should
      contain each unique allele for the genotype once.  If an array
      element contains a vector of length 1 containing only the symbol
      that is in the \code{Missing} slot, this indicates missing data
      for that sample and locus.}
    \item{\code{Description}:}{Object of class \code{"character"}.  This
    stores a description of the dataset for the user's convenience.}
    \item{\code{Missing}:}{Object of class \code{"ANY"}.  A symbol to be
    used to indicate missing data in the \code{Genotypes} slot.  This is
  the integer \code{-9} by default.}
    \item{\code{Usatnts}:}{Object of class \code{"integer"}.  A vector,
      named by loci.  Each element indicates the repeat type of the
      locus.  \code{2} indicates dinucleotide repeats, \code{3}
      indicates trinucleotide repeats, and so on.  If the alleles stored
    in the \code{Genotypes} slot for a given locus are already written
    in terms of repeat number, the \code{Usatnts} value for that locus
    should be \code{1}.  In other words, all alleles for a locus can be
    divided by the number in \code{Usatnts} to give alleles expressed in
  terms of relative repeat number.}
    \item{\code{Ploidies}:}{Object of class \code{"integer"}.  A vector,
    named by samples.  This stores the ploidy of each sample.  \code{NA}
  indicates unknown ploidy.  See \code{\link{Ploidies<-}} and
  \code{\link{estimatePloidy}} for ways to fill this slot.}
    \item{\code{PopInfo}:}{Object of class \code{"integer"}.  A vector,
      named by samples, containing the population identity of each sample. }
    \item{\code{PopNames}:}{Object of class \code{"character"}.  A
      vector containing names for all populations.  The position of a
      population name in the vector indicates the integer used to
      represent that population in \code{PopInfo}. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gendata}"}, directly.
}
\section{Methods}{
  For more information on any of these methods, see the help files of
  their respective generic functions.
  \describe{
    \item{deleteLoci}{\code{signature(object = "genambig")}: Removes
      columns in the array in the \code{Genotypes} slot corresponding to
    the locus names supplied, then passes the arguments to the method
    for \code{gendata}.}
    \item{deleteSamples}{\code{signature(object = "genambig")}: Removes
      rows in the array in the \code{Genotypes} slot corresponding to
      the sample names supplied, then passes the arguments to the method
    for \code{gendata}.}
    \item{editGenotypes}{\code{signature(object = "genambig")}: Each
      vector in the \code{Genotypes} slot is placed into the row of a
      data frame, along with the sample and locus name for this vector.
    The data frame is then opened in the Data Editor so that the user
    can make changes.  When the Data Editor window is closed, vectors
    are extracted back out of the data frame and written to the
    \code{Genotypes} slot.}
    \item{estimatePloidy}{\code{signature(object = "genambig")}:
      Calculates the length of each genotype vector (excluding those
      with the missing data symbol), and creates a data frame showing
      the maximum and mean number of alleles per locus for each sample.
    This data frame is then opened in the Data Editor, where the user
    may edit ploidy levels.  Once the Data Editor is closed, the
    \code{genambig} object is returned with the new values written to
    the \code{Ploidies} slot.}
    \item{Genotype}{\code{signature(object = "genambig")}: Retrieves a
      single genotype vector, as specified by \code{sample} and
      \code{locus} arguments. }
    \item{Genotype<-}{\code{signature(object = "genambig")}: Replaces a
      single genotype vector. }
    \item{Genotypes}{\code{signature(object = "genambig")}: Retrieves a
      two-dimensional list of genotype vectors. }
    \item{Genotypes<-}{\code{signature(object = "genambig")}: Replaces a
    one- or two-dimensional list of genotype vectors.}
    \item{initialize}{\code{signature(.Object = "genambig")}: When
      \code{new} is called to create a new \code{genambig} object, the
      \code{initialize} method sets up a two dimensional list in the
      \code{Genotypes} slot indexed by sample and locus, and fills this
      list with the missing data symbol.  The \code{initialize} method
      for \code{gendata} is then called. }
    \item{isMissing}{\code{signature(object = "genambig")}: Given a set of
    samples and loci, each position in the array in the \code{Genotypes}
  slot is checked to see if it matches the missing data value.  A single
Boolean value or an array of Boolean values is returned.}
    \item{Loci<-}{\code{signature(object = "genambig")}: For changing
      the names of loci.  The names are changed in the second dimension
      of the array in the \code{Genotypes} slot, and then the
      \code{Loci<-} method for \code{gendata} is called. }
    \item{Missing<-}{\code{signature(object = "genambig")}: For changing
    the missing data symbol.  All elements of the \code{Genotypes} array
  that match the current missing data symbol are changed to the new
  missing data symbol.  The \code{Missing<-} method for \code{gendata}
  is then called.}
    \item{Samples<-}{\code{signature(object = "genambig")}: For changing
    the names of samples.  The names are changed in the first dimension
    of the array in the \code{Genotypes} slot, and then the
    \code{Samples<-} method for \code{gendata} is called.}
    \item{summary}{\code{signature(object = "genambig")}: Prints the
      dataset description (\code{Description} slot) to the console as
      well as the number of missing genotypes, then calls the
      \code{summary} method for \code{gendata}. }
    \item{show}{\code{signature(object = "genambig")}: Prints the
      data to the console, formatted to make it more legible. The
      genotype for each locus is shown as the size of each allele,
      separated by a '/'. The SSR motif length ('Usatnts'), ploidies,
      population names, and population membership ('PopInfo') are
      displayed if they exist.}
    \item{viewGenotypes}{\code{signature(object = "genambig")}: Prints a
    tab-delimited table of samples, loci, and genotype vectors to the console.}
    \item{"["}{\code{signature(x = "genambig", i = "ANY", j = "ANY")}:
      For subscipting \code{genambig} objects.  Should be of the form
      \code{mygenambig[mysamples, myloci]}.  Returns a \code{genambig}
      object.  The \code{Genotypes} slot
      is replaced by one containing only samples \code{i} and loci
      \code{j}.  Likewise, the \code{PopInfo} and \code{Ploidies} slots
      are truncated to contain only samples \code{i}, and the
      \code{Usatnts} slot is truncated to contain only loci
      \code{j}.  Other slots are left unaltered.}
    \item{merge}{\code{signature(x = "genambig", y = "genambig")}:
      Merges two genotypes objects together.  See
      \code{\link{merge,genambig,genambig-method}}.}
	 }
}

\author{
Lindsay V. Clark, Tyler W. Smith
}


\seealso{
\code{\linkS4class{gendata}}, \code{\link[=Samples]{Accessors}},
\code{\link{merge,genambig,genambig-method}}
}
\examples{
# display class definition
showClass("genambig")

# create a genambig object
mygen <- new("genambig", samples=c("a", "b", "c", "d"),
             loci=c("L1", "L2", "L3"))
# add some genotypes
Genotypes(mygen)[,"L1"] <- list(c(133, 139, 145), c(142, 154),
                                c(130, 142, 148), Missing(mygen))
Genotypes(mygen, loci="L2") <- list(c(105, 109, 113), c(111, 117),
                                    c(103, 115), c(105, 109, 113))
Genotypes(mygen, loci="L3") <- list(c(254, 258), Missing(mygen),
                                    c(246, 250, 262), c(250, 258))

# see a summary of the object
summary(mygen)
# display some of the genotypes
viewGenotypes(mygen[c("a", "b", "c"),])
}
\keyword{classes}
