\name{estimatePloidy}
\alias{estimatePloidy}
\title{
Estimate Ploidies Based on Allele Counts
}
\description{
\code{estimatePloidy} calculates the maximum and mean number of unique
alleles for each sample across a given set of loci.  These values are
presented in a data editor, along with other pertinent information, so
that the user can then edit the ploidy values for the object.
}
\usage{
estimatePloidy(object, extrainfo, samples = Samples(object), 
               loci = Loci(object))
}
\arguments{
  \item{object}{
The object containing genotype data, and to which ploidies will be written.
}
  \item{extrainfo}{
A named or unnamed vector or data frame containing extra information
(such as morphological or flow cytometry data) to
display in the data editor, to assist with making decisions about
ploidy.  If unnamed, the vector (or the rows of the data frame) is
assumed to be in the same order as \code{samples}.  An array can also be
given as an argument here, and will be coerced to a data frame.
}
  \item{samples}{
A numeric or character vector indicating a subset of samples to evaluate.
}
  \item{loci}{
A numeric or character vector indicating a subset of loci to use in the
calculation of mean and maximum allele number.
}
}
\details{
\code{estimatePloidy} is a generic function with methods written for
the \code{genambig} and \code{genbinary} classes.

If the \code{Ploidies} slot of \code{object} is not already a
\code{"ploidysample"} object, the function will first convert the
\code{Ploidies} slot to this format, deleting any data that is currently
there.  (Ploidies must be indexed by sample and not by locus.)  If
ploidies were already in the \code{"ploidysample"} format, any ploidy
data already in the object is retained and put into the table (see
below).

Population identities are displayed in the table only if more than one
population identity is found in the dataset.  Likewise, the current
ploidies of the dataset are only displayed if there is more than one
ploidy level already found in \code{Ploidies(object)}.

Missing genotypes are ignored; maximum and mean allele counts are only
calculted across genotypes that are not missing.  If all genotypes for a
given sample are missing, \code{NA} is displayed in the corresponding
cells in the data editor.

The default values for \code{new.ploidy} are the maximum number of
alleles per locus for each sample.
}
\value{
\code{object} is returned, with \code{Ploidies(object)} now equal to the
values set in the \code{new.ploidy} column of the data editor.
}

\author{
Lindsay V. Clark
}


\seealso{
\code{\linkS4class{genambig}}, \code{\linkS4class{genbinary}},
\code{\link{Ploidies}}
}
\examples{
if(interactive()){ #this line included for automated checking on CRAN

# create a dataset for this example
mygen <- new("genambig", samples=c("a", "b", "c"),
             loci=c("loc1", "loc2"))
Genotypes(mygen, loci="loc1") <- list(c(122, 126, 128), c(124, 130),
                                      c(120, 122, 124))
Genotypes(mygen, loci="loc2") <- list(c(140, 148), c(144, 150), Missing(mygen))

# estimate the ploidies
mygen <- estimatePloidy(mygen)

# view the ploidies
Ploidies(mygen)

}

}

\keyword{ arith }
\keyword{ methods }
