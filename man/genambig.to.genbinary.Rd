\name{genambig.to.genbinary}
\alias{genambig.to.genbinary}
\alias{genbinary.to.genambig}
\title{
Convert Between Genotype Object Classes
}
\description{
These functions convert back and forth between the \code{genambig} and
\code{genbinary} classes.
}
\usage{
genambig.to.genbinary(object, samples = Samples(object),
                      loci = Loci(object))

genbinary.to.genambig(object, samples = Samples(object),
                      loci = Loci(object))
}

\arguments{
  \item{object}{
The object containing the genetic dataset.  A \code{genambig} object for
  \code{genambig.to.genbinary}, or a \code{genbinary} object for
  \code{genbinary.to.genambig}.
}
  \item{samples}{
An optional character vector indicating samples to include in the new object.
}
  \item{loci}{
An optional character vector indicating loci to include in the new object.
}
}
\details{
The slots \code{Description}, \code{Ploidies}, \code{Usatnts},
  \code{PopNames}, and \code{PopInfo} are transferred as-is from the old
  object to the new.  The value in the
  \code{Genotypes} slot is converted from one format to the other, with
  preservation of allele names.
}
\value{
For \code{genambig.to.genbinary}: a \code{genbinary} object containing
  all of the data from \code{object}.  \code{Missing}, \code{Present},
  and \code{Absent} are set at their default values.

  For \code{genbinary.to.genambig}: a \code{genambig} object containing
  all of the data from \code{object}.  \code{Missing} is at the default
  value.
}

\author{
Lindsay V. Clark
}


\seealso{
\code{\linkS4class{genambig}}, \code{\linkS4class{genbinary}}
}
\examples{
# set up a genambig object for this example
mygen <- new("genambig", samples = c("A", "B", "C", "D"),
             loci = c("locJ", "locK"))
PopNames(mygen) <- c("PopQ", "PopR")
PopInfo(mygen) <- c(1,1,2,2)
Usatnts(mygen) <- c(2,2)
Genotypes(mygen, loci="locJ") <- list(c(178, 184, 186), c(174,186),
                                      c(182, 188, 190),
                                      c(182, 184, 188))
Genotypes(mygen, loci="locK") <- list(c(133, 135, 141),
                                      c(131, 135, 137, 143),
                                      Missing(mygen), c(133, 137))

# convert it to a genbinary object
mygenB <- genambig.to.genbinary(mygen)

# check the results
viewGenotypes(mygenB)
viewGenotypes(mygen)
PopInfo(mygenB)

# convert back to a genambig object
mygenA <- genbinary.to.genambig(mygenB)
viewGenotypes(mygenA)

# note: identical(mygen, mygenA) returns FALSE, because the alleles
# origninally input are not stored as integers, while the alleles
# produced by genbinary.to.genambig are integers.

  }


\keyword{ manip }

