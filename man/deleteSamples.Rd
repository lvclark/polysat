\name{deleteSamples}
\alias{deleteSamples}
\alias{deleteLoci}
\title{
Remove Samples or Loci from an Object
}
\description{
These functions remove samples or loci from all relevant slots of an object.
}
\usage{
deleteSamples(object, samples)
deleteLoci(object, loci)
}
\arguments{
  \item{object}{
An object containing the dataset of interest.  Generally an object of
some subclass of \code{gendata}.
}
  \item{samples}{
A numerical or character vector of samples to be removed.
}
\item{loci}{
  A numerical or character vector of loci to be removed.
  }
}
\details{
These are generic functions with methods for \code{genambig},
\code{genbinary}, and
\code{gendata} objects.  The methods for the subclasses remove samples
or loci
from the \code{@Genotypes} slot, then pass the object to the method for
\code{gendata}, which removes samples or loci from the \code{@PopInfo},
\code{@Ploidies}, and/or \code{@Usatnts} slots, as appropriate.  The
\code{@PopNames} slot is left untouched even if an entire population is
deleted, in order to preserve the connection between the numbers in
\code{@PopInfo} and the names in \code{@PopNames}.

If your intent is to experiment with excluding samples or loci, it may
be a better idea to create character vectors of samples and loci that
you want to use and then use these vectors as the \code{samples} and
\code{loci} arguments for analysis or export functions.
}
\value{
An object identical to \code{object}, but with the specified samples or
loci removed.
}

\author{
Lindsay V. Clark
}
\note{
These functions are somewhat redundant with the subscripting function
\code{"["}, which also works for all \code{gendata} objects.  However,
they may be more convenient depending on whether the user prefers to
specify the samples and loci to use or to
exclude.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{Samples}}, \code{\link{Loci}},
  \code{\link{merge,gendata,gendata-method}}
}
\examples{
# set up genambig object
mygen <- new("genambig", samples = c("ind1", "ind2", "ind3", "ind4"),
             loci = c("locA", "locB", "locC", "locD"))

# delete a sample
Samples(mygen)
mygen <- deleteSamples(mygen, "ind1")
Samples(mygen)

# delete some loci
Loci(mygen)
mygen <- deleteLoci(mygen, c("locB", "locC"))
Loci(mygen)
}
\keyword{ manip }

