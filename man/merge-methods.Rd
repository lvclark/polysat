\name{merge-methods}
\docType{methods}
\alias{merge}
\alias{merge-methods}
\alias{merge,genambig,genambig-method}
\alias{merge,gendata,gendata-method}
\alias{merge,genbinary,genbinary-method}
\title{ Merge Two Genotype Objects into One }
\description{
The generic function \code{merge} has methods defined in polysat to
merge two genotype objects of the same class.  Each method has optional
\code{samples} and \code{loci} arguments for specifying subsets of
samples and loci to be included in the merged object.  Each method also
has an optional \code{overwrite} argument to specify which of the two
objects should not be used in the case of conflicting data.
}
\usage{
merge(x, y, ...)
}
\arguments{
  \item{x}{One of the objects to be merged.  For the methods defined for
    \pkg{polysat} this should be of class \code{"gendata"} or one of its
    subclasses.
  }
  \item{y}{
    The other object to be merged.  Should be of the same class as
    \code{x}.  \code{y@Ploidies} must also be of the same class as
    \code{x@Ploidies}.
  }
  \item{...}{
    Additional arguments specific to the method.
    }
  }
  \section{Methods}{

    The methods for \code{merge} in polysat have four additional
    arguments: \code{objectm, samples, loci, overwrite}.

    The \code{samples} and \code{loci} arguments can specify, using
    character vectors, a subset of the samples and loci found in
    \code{x} and \code{y} to write to the object that is returned.

    If \code{overwrite = "x"}, data from the second object will be used
    wherever there is contradicting data.  Likewise if \code{overwrite =
    "y"}, data from the first object will be used wherever there is
    contradicting data.  If no \code{overwrite} argument is given, then
    any contradicting data between the two objects will produce an error
    indicating where the contradicting data were found.

    The \code{objectm} argument is primarily for internal use (most
    users will not need it).  If this argument is not
    provided, a new genotype object is created and data from \code{x}
    and \code{y} are written to it.  If \code{objectm} is
    provided, this is the object to which data will be written, and the
    object that will be returned.

  \describe{
\item{\code{signature(x = "genambig", y = "genambig")}}{
This method merges the genotype data from \code{x} and \code{y}.  If the
missing data symbols differ between the objects, \code{overwrite} is
used to determine which missing data symbol to use, and all missing data
symbols in the overwritten object are converted.  If \code{overwrite} is
not provided and the missing data symbols differ between the objects, an
error will be given.  The genotypes are then filled in.  If certain
sample*locus combinations do not exist in either object (\code{x} and
\code{y} have different samples as well as different loci), missing data
symbols are left in these positions.  Again, for genotypes,
\code{overwrite} determines which object to preferentially use for data
and whether to give an error if there is a disagreement.

The \code{merge} method for \code{gendata} is then called.
}

\item{\code{signature(x = "genbinary", y = "genbinary")}}{
This method also merges genotype data for \code{x} and \code{y}, then
calls the method for \code{gendata}.  \code{Missing}, \code{Present},
and \code{Absent} are checked for consistency between objects similarly
to what happens with \code{Missing} in the \code{genambig} method.
Genotypes are then written to the merged object, and consistency between
genotypes is checked.
}

\item{\code{signature(x = "gendata", y = "gendata")}}{
This method merges data about ploidy, repeat length, and population
identity, as well as writing one or both dataset descriptions to the
merged object.

The same population numbers can have different meanings in
\code{PopInfo(x)} and \code{PopInfo(y)}.  The unique \code{PopNames} are
used instead to determine population identity, and the \code{PopInfo}
numbers are changed if necessary.  Therefore, it is important for
identical populations to be named the same way in both objects, but not
important for identical populations to have the same number in both objects.
}
}}
\keyword{methods}
\keyword{ manip }
