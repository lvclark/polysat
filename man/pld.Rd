\name{pld}
\alias{pld}
\alias{pld<-}
\alias{plCollapse}
\title{
Accessor, Replacement, and Manipulation Functions for \code{"ploidysuper"} Objects
}
\description{
\code{pld} accesses and replaces the \code{pld} slot of objects of 
\code{"\linkS4class{ploidysuper}"} subclasses.  \code{plCollapse} tests
whether an object of one of these classes can be converted to an object
of a simpler one of these classes, and optionally returns the converted
object.  These are generic functions with methods for the subclasses of
\code{"ploidysuper"}.  These functions are primarily for internal use.
}
\usage{
pld(object, samples, loci)
pld(object) <- value
plCollapse(object, na.rm, returnvalue)
}
\arguments{
  \item{object}{
A \code{"ploidysuper"} object.
}
  \item{samples}{
An optional character or numeric vector indexing the samples for which
to return ploidy values.
}
\item{loci}{
  An optional character or numeric vector indexing the loci for which
to return ploidy values.
}
\item{value}{
  A numeric vector or matrix that can be coerced to integers.  These
  represent the ploidies to store in the \code{object@pld} slot.
}
\item{na.rm}{Boolean.  If \code{TRUE}, \code{NA} values are ignored when
  testing to see if the ploidy format can be simplified.  If the sample,
  locus, or entire dataset all has one ploidy aside from \code{NA}, the
  \code{NA} values will be overwritten by that ploidy when simplifying
  the ploidy format.  If \code{FALSE}, \code{NA} is treated as a unique
  ploidy.
}
\item{returnvalue}{Boolean.  If \code{TRUE}, a \code{"ploidysuper"}
  object will be returned if the ploidy format can be simplified, and
  \code{FALSE} will be returned if it cannot be simplified.  If
  \code{FALSE}, only \code{TRUE} or \code{FALSE} will be returned to
  indicate if the ploidy format can be simplified or not.
  }
}

\value{
\code{pld} returns the vector or matrix containing the ploidy values.
This is the contents of \code{object@pld}.

\code{plCollapse} either returns a Boolean value indicating whether the
ploidy can be changed to a simpler format, or a new \code{"ploidysuper"}
object with all of the ploidy data of \code{object} put into a simpler
format.  If \code{object} is a \code{"ploidymatrix"} object, a
\code{"ploidysample"}, \code{"ploidylocus"}, or \code{"ploidyone"}
object can be returned depending on how many unique ploidy values there
are and how they are distributed.  If \code{object} is a
\code{"ploidysample"} or \code{"ploidylocus"} object, a
\code{"ploidyone"} object can be returned.
}

\author{
Lindsay V. Clark
}


\seealso{
\code{\link{reformatPloidies}}, \code{\link{Ploidies}}
}
\examples{
test <- new("ploidymatrix", samples=c("a","b","c"),
            loci=c("l1","l2","l3"))
pld(test)       # view the ploidies
pld(test) <- 2  # make it diploid at all samples and loci
pld(test)["a",] <- c(2,4,4) # change the ploidies for sample a
pld(test, samples=c("a","b")) # view ploidies at a subset of samples

# test to see if the ploidies can be simplified
p <- plCollapse(test, na.rm=FALSE, returnvalue=TRUE)
p
# now change a ploidy and repeat the test
pld(test)["a","l1"] <- 4
p <- plCollapse(test, na.rm=FALSE, returnvalue=TRUE)
p
# change something else and collapse it further
pld(p)["a"] <- 2
p2 <- plCollapse(p, na.rm=FALSE, returnvalue=TRUE)
p2

# if na.rm=FALSE, NA values are not ignored:
pld(test)["a","l1"] <- NA
pld(test)
plCollapse(test, na.rm=FALSE, returnvalue=TRUE)
# NA values are ignored with na.rm=TRUE
plCollapse(test, na.rm=TRUE, returnvalue=TRUE)
}
\keyword{ methods }
\keyword{ manip }% __ONLY ONE__ keyword per line
