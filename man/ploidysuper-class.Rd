\name{ploidysuper-class}
\Rdversion{1.1}
\docType{class}
\alias{ploidysuper-class}
\alias{ploidymatrix-class}
\alias{ploidysample-class}
\alias{ploidylocus-class}
\alias{ploidyone-class}
\alias{pld,ploidymatrix-method}
\alias{pld,ploidysample-method}
\alias{pld,ploidylocus-method}
\alias{pld,ploidyone-method}
\alias{pld<-,ploidymatrix-method}
\alias{pld<-,ploidysample-method}
\alias{pld<-,ploidylocus-method}
\alias{pld<-,ploidyone-method}
\alias{plCollapse,ploidymatrix,logical,logical-method}
\alias{plCollapse,ploidysample,logical,logical-method}
\alias{plCollapse,ploidylocus,logical,logical-method}
\alias{plCollapse,ploidyone,logical,logical-method}

\title{Class \code{"ploidysuper"} and Subclasses}
\description{
These classes contain ploidy data indexed by sample, locus, both, or
neither.  They are intended to go in the \code{\link{Ploidies}} slot of
\code{"\linkS4class{gendata}"} objects.
}
\section{Objects from the Class}{\code{"ploidysuper"} is a virtual
  class: No objects may be created from it.

Objects of the subclasses \code{"ploidymatrix"}, \code{"ploidysample"},
\code{"ploidylocus"}, and \code{"ploidyone"} can be created with the
call \code{new(ploidyclass, samples, loci, ...)}, where
\code{ploidyclass} is a character string of one of the class names, and
\code{samples} and \code{loci} are character vectors naming samples and
loci, respectively.  The latter two arguments are optional depending on
the class (whether ploidies are indexed by sample and/or locus).  The
typical user will not have to create an object in this way, because other
functions in \pkg{polysat} will do it for you.}
\section{Slots}{
  \describe{
    \item{\code{pld}:}{The only slot for objects of these classes.  For
      \code{"ploidymatrix"}, this is a matrix of integers, indexed in
      the first dimension by sample name and in the second dimension by
      locus name.  Each element represents the ploidy at a given sample
      and locus.  For \code{"ploidysample"} and \code{"ploidylocus"},
      the slot is an integer vector, named by sample or locus and
      indicating the ploidy at each sample or locus, respectively.  For
      \code{"ploidyone"}, the slot contains a single integer
      representing the ploidy for the entire dataset.}
    }
}
\section{Methods}{
  \describe{
\item{pld}{\code{signature(object= "ploidymatrix")}: Returns the
  contents of \code{object@pld}: a matrix of ploidies indexed by sample
  and locus.  The \code{samples} and \code{loci} arguments can be used,
  optionally, to only return a subset of ploidies.}
}}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{reformatPloidies}}, \code{\link{pld}},
\code{\link{plCollapse}}, \code{\linkS4class{gendata}}
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("ploidysuper")
}
\keyword{classes}
