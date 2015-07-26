\name{gendata-class}
\Rdversion{1.1}
\docType{class}
\alias{gendata-class}
\alias{deleteLoci,gendata-method}
\alias{deleteSamples,gendata-method}
\alias{Description,gendata-method}
\alias{Description<-,gendata-method}
\alias{initialize,gendata-method}
\alias{Loci,gendata,missing,missing-method}
\alias{Loci,gendata,numeric,missing-method}
\alias{Loci,gendata,missing,numeric-method}
\alias{Loci,gendata,numeric,numeric-method}
\alias{Loci<-,gendata-method}
\alias{Missing,gendata-method}
\alias{Missing<-,gendata-method}
\alias{Ploidies,gendata-method}
\alias{Ploidies<-,gendata-method}
\alias{PopInfo,gendata-method}
\alias{PopInfo<-,gendata-method}
\alias{PopNames,gendata-method}
\alias{PopNames<-,gendata-method}
\alias{PopNum,gendata,character-method}
\alias{PopNum<-,gendata,character-method}
\alias{Samples,gendata,character,missing-method}
\alias{Samples,gendata,character,numeric-method}
\alias{Samples,gendata,missing,missing-method}
\alias{Samples,gendata,missing,numeric-method}
\alias{Samples,gendata,numeric,missing-method}
\alias{Samples,gendata,numeric,numeric-method}
\alias{Samples<-,gendata-method}
\alias{summary,gendata-method}
\alias{Usatnts,gendata-method}
\alias{Usatnts<-,gendata-method}
\alias{[,gendata-method}

\title{Class "gendata"}
\description{
This is a superclass for other classes that contain population
genetic datasets.  It has slots for population information, ploidy,
microsatellite repeat lengths, and a missing data symbol, but does not
have a slot to store genotypes.  Sample and locus names are stored as
the \code{names} of vectors in the slots.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gendata", samples, loci, ...)}.
The missing data symbol will be set to \code{-9} by default.  The
default initial value for \code{PopNames} is a
character vector of length 0, and for \code{Description} is the string
\code{"Insert dataset description here"}.  The default initial value for
the \code{Ploidies} slot is a \code{"\linkS4class{ploidymatrix}"}
object, containing a matrix filled with \code{NA} and named by samples
in the first dimension and loci in the second dimension.  For
other slots, vectors filled with \code{NA} will be generated and will be
named by samples (for \code{PopInfo}) or loci (for
\code{Usatnts}).  The slots can then be edited
using the methods described below.

Note that in most cases you will want to instead create an
object from one of \code{gendata}'s subclasses, such as \code{genambig}.
}
\section{Slots}{
  \describe{
    \item{\code{Description}:}{Object of class \code{"character"}.  One
      or more character strings to name or describe the dataset. }
    \item{\code{Missing}:}{Object of class \code{"ANY"}.  A value to
      indicate missing data in the genotypes of the dataset.  \code{-9}
      by default. }
    \item{\code{Usatnts}:}{Object of class \code{"integer"}.  This
      vector must be named by locus names.  Each element should be the
      length of the microsatellite repeat for that locus, given in nucleotides.
     For example, \code{2} would indicate a locus with dinucleotide
     repeats, and \code{3} would indicate a locus with trinucleotide
     repeats.  \code{1} should be used for mononucleotide repeats OR if
     alleles for that locus are already expressed in terms of repeat
     number rather than nucleotides.  To put it another way, if you
     divided the number used to represent an allele by the corresponding
   number in \code{Usatnts} (and rounded if necessary), the result would
 be the number of repeats (plus some additional length for flanking regions). }
    \item{\code{Ploidies}:}{Object of class \code{"ploidysuper"}.  This
      object will contain a \code{pld} slot that is a matrix named by
      samples and loci, a vector named by samples or loci, or a single
      value, depending on the subclass.  Each element is an integer that
      represents ploidy.  \code{NA} indicates unknown ploidy. }
    \item{\code{PopInfo}:}{Object of class \code{"integer"}.  This
      vector also must be named by sample names.  Each element
      represents the number of the population to which each sample belongs. }
    \item{\code{PopNames}:}{Object of class \code{"character"}.  An
      unnamed vector containing the name of each population.  If a
      number from \code{PopInfo} is used to index \code{PopNames}, it
      should find the correct population name.  For example, if the
      first element of \code{PopNames} is \code{"ABC"}, then any samples
    with \code{1} as their \code{PopInfo} value belong to population
    \code{"ABC"}.}
  }
}
\section{Methods}{
  \describe{
    \item{deleteLoci}{\code{signature(object = "gendata")}: Permanently
      remove loci from the dataset.  This removes elements from \code{Usatnts}. }
    \item{deleteSamples}{\code{signature(object = "gendata")}:
      Permanently remove samples from the dataset.  This removes
      elements from \code{PopInfo} and \code{Ploidies}. }
    \item{Description}{\code{signature(object = "gendata")}: Returns the
    character vector in the \code{Description} slot.}
    \item{Description<-}{\code{signature(object = "gendata")}: Assigns a
    new value to the character vector in the \code{Description} slot.}
    \item{initialize}{\code{signature(.Object = "gendata")}: This is
      called when the \code{new("gendata")} function is used.  A new
      \code{gendata} object is created with sample and locus names used
      to index the appropriate slots. }
    \item{Loci}{\code{signature(object = "gendata", usatnts =
	"missing", ploidies="missing")}: Returns a character vector
      containing all locus
      names for the object.  The method accomplishes this by returning
      \code{names(object@Usatnts)}. }
    \item{Loci}{\code{signature(object = "gendata", usatnts =
	"numeric", ploidies = "missing")}: Returns a character vector
      of all loci for a given
      set of repeat lengths.  For example, if \code{usatnts = 2} all
      loci with dinucleotide repeats will be returned.}
    \item{Loci}{\code{signature(object= "gendata", usatnts = "missing",
    ploidies = "numeric")}: Returns a character vector of all loci for a
  given set of ploidies. Only works if \code{object@Ploidies} is a
  \code{"\linkS4class{ploidylocus}"} object.}
\item{Loci}{\code{signature(object = "gendata", usatnts = "numeric",
ploidies = "numeric")}: Returns a character vector of all loci that have
one of the indicated repeat types and one of the indicated ploidies.
Only works if \code{object@Ploidies} is a \code{"ploidylocus"} object.}
    \item{Loci<-}{\code{signature(object = "gendata")}: Assigns new
      names to loci in the dataset (changes \code{names(object@Usants)}.
     Should not be used for adding or removing loci.}
    \item{Missing}{\code{signature(object = "gendata")}: Returns the
      missing data symbol from \code{object@Missing}. }
    \item{Missing<-}{\code{signature(object = "gendata")}: Assigns a new
    value to \code{object@Missing} (changes the missing data symbol).}
  \item{Ploidies}{\code{signature(object = "gendata", samples = "ANY",
    loci = "ANY")}: Returns the
      ploidies in the dataset (\code{object@Ploidies}), indexed by
      sample and locus if applicable.. }
    \item{Ploidies<-}{\code{signature(object = "gendata")}: Assigns new
      values to ploidies of samples in the dataset. The assigned values
      are coerced to integers by the method.  Names in the assigned
      vector or matrix are ignored; sample and/or locus names already
      present in the \code{gendata} object are used instead.}
    \item{PopInfo}{\code{signature(object = "gendata")}: Returns the
      population numbers of samples in the dataset  (\code{object@PopInfo}). }
    \item{PopInfo<-}{\code{signature(object = "gendata")}: Assigns new
      population numbers to samples in the dataset.  The assigned values
      are coerced to integers by the method.  Names in the assigned
      vector are ignored; sample names already present in the
      \code{gendata} object are used instead. }
    \item{PopNames}{\code{signature(object = "gendata")}: Returns a
      character vector of population names (\code{object@PopNames}). }
    \item{PopNames<-}{\code{signature(object = "gendata")}: Assigns new
      names to populations. }
    \item{PopNum}{\code{signature(object = "gendata",
	popname="character")}: Returns the number corresponding to a
      population name.}
    \item{PopNum<-}{\code{signature(object = "gendata", popname =
	"character")}: Changes the population number for a given
      population name, merging it with an existing population of that
      number if applicable.}
    \item{Samples}{\code{signature(object = "gendata", populations =
	"character", ploidies = "missing")}: Returns all sample names
      for a given set of population names. }
    \item{Samples}{\code{signature(object = "gendata", populations =
	"character", ploidies = "numeric")}: Returns all sample names
      for a given set of population names and ploidies.  Only samples
      that fit both criteria will be returned. }
    \item{Samples}{\code{signature(object = "gendata", populations =
	"missing", ploidies = "missing")}: Returns all sample names. }
    \item{Samples}{\code{signature(object = "gendata", populations =
	"missing", ploidies = "numeric")}: Returns all sample names for
      a given set of ploidies. Only works if \code{object@Ploidies} is a
    \code{"\linkS4class{ploidysample}"} object.}
    \item{Samples}{\code{signature(object = "gendata", populations =
	"numeric", ploidies = "missing")}: Returns all sample names for
      a given set of population numbers. }
    \item{Samples}{\code{signature(object = "gendata", populations =
	"numeric", ploidies = "numeric")}: Returns all sample names for
      a given set of population numbers and ploidies.  Only samples that
    fit both criteria will be returned.  Only works if
    \code{object@Ploidies} is a \code{"ploidysample"} object.}
    \item{Samples<-}{\code{signature(object = "gendata")}: Assigns new
      names to samples.  This edits both \code{names(object@PopInfo)} and
      \code{names(object@Ploidies)}.  It should not be used for adding
      or removing samples from the dataset.}
    \item{summary}{\code{signature(object = "gendata")}: Prints some
      informaton to the console, including the numbers of samples, loci,
    and populations, the ploidies present, and the types of
    microsatellite repeats present.}
    \item{Usatnts}{\code{signature(object = "gendata")}: Returns
      microsatellite repeat lengths for loci in the dataset
      (\code{object@Usatnts}). }
    \item{Usatnts<-}{\code{signature(object = "gendata")}: Assigns new
      values to microsatellite repeat lengths of loci
      (\code{object@Usatnts}).  The assigned values
      are coerced to integers by the method.  Names in the assigned
      vector are ignored; locus names already present in the
      \code{gendata} object are used instead. }
    \item{"["}{\code{signature(x = "gendata", i = "ANY", j = "ANY")}:
      Subscripts the data by a subset of samples and/or loci.  Should be
      used in the format \code{mygendata[mysamples, myloci]}.  Returns a
      \code{gendata} object with \code{PopInfo}, \code{Ploidies}, and
      \code{Usatnts} truncated to only contain the samples and loci
      listed in \code{i} and \code{j}, respectively.
      \code{Description}, \code{Missing}, and \code{PopNames} are left
      unaltered.}
    \item{merge}{\code{signature(x = "gendata", y = "gendata")}: Merges
      two genotype objects.  See \code{\link{merge,gendata,gendata-method}}.}
	 }
}

\author{
Lindsay V. Clark
}

\seealso{
  \code{\linkS4class{genambig}}, \code{\linkS4class{genbinary}},
  \code{\link[=Samples]{Accessors}}
}
\examples{
# show class definition
showClass("gendata")

# create an object of the class gendata
# (in reality you would want to create an object belonging to one of the
# subclasses, but the procedure is the same)
mygen <- new("gendata", samples = c("a", "b", "c"),
             loci = c("loc1", "loc2"))
Description(mygen) <- "An example for the documentation"
Usatnts(mygen) <- c(2,3)
PopNames(mygen) <- c("PopV", "PopX")
PopInfo(mygen) <- c(2,1,2)
Ploidies(mygen) <- c(2,2,4,2,2,2)

# view a summary of the object
summary(mygen)
}
\keyword{classes}
