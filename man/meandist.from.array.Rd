\name{meandist.from.array}
\alias{meandist.from.array}
\alias{find.na.dist}
\alias{find.na.dist.not.missing}
\title{
Tools for Working With Pairwise Distance Arrays
}
\description{
\code{meandist.from.array} produces a mean distance matrix from an array of
pairwise distances by locus, such as that produced by
\code{meandistance.matrix} when \code{all.distances=TRUE}.  \code{find.na.dist}
finds missing distances in such an array, and
\code{find.na.dist.not.missing} finds missing distances that aren't the
result of missing genotypes.
}
\usage{
meandist.from.array(distarray, samples = dimnames(distarray)[[2]],
loci = dimnames(distarray)[[1]])

find.na.dist(distarray, samples = dimnames(distarray)[[2]],
loci = dimnames(distarray)[[1]])

find.na.dist.not.missing(object, distarray,
samples = dimnames(distarray)[[2]], loci = dimnames(distarray)[[1]])
}
\arguments{
  \item{distarray}{
A three-dimensional array of pairwise distances between samples, by
locus.  Loci are represented in the first dimension, and samples are
represented in the second and third dimensions.  Dimensions are named
accordingly.  Such an array is the first element of the list produced by
\code{meandistance.matrix} if \code{all.distances=TRUE}.
}
  \item{samples}{
Character vector.  Samples to analyze.
}
  \item{loci}{
Character vector.  Loci to analyze.
}
\item{object}{
  A \code{genambig} object.  Typically the
  genotype object that was used to produce \code{distarray}.
}
}
\details{
\code{find.na.dist.not.missing} is primarily intended to locate distances
that were not calculated by \code{Bruvo.distance} because both genotypes
had too many alleles (more than \code{maxl}).  The user may wish to
estimate these distances manually and fill them into the array, then
recalculate the mean matrix using \code{meandist.from.array}.
}
\value{
\code{meandist.from.array} returns a matrix, with both rows and columns
named by samples, of distances averaged across loci.

\code{find.na.dist} and \code{find.na.dist.not.missing} both return data
frames with three columns: Locus, Sample1, and Sample2.  Each row
represents the index in the array of an element containing NA.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{meandistance.matrix}}, \code{\link{Bruvo.distance}},
\code{\link{find.missing.gen}}
}
\examples{
# set up the genotype data
samples <- paste("ind", 1:4, sep="")
samples
loci <- paste("loc", 1:3, sep="")
loci
testgen <- new("genambig", samples=samples, loci=loci)
Genotypes(testgen, loci="loc1") <- list(c(-9), c(102,104),
                                        c(100,106,108,110,114),
                                        c(102,104,106,110,112))
Genotypes(testgen, loci="loc2") <- list(c(77,79,83), c(79,85), c(-9),
                                        c(83,85,87,91))
Genotypes(testgen, loci="loc3") <- list(c(122,128), c(124,126,128,132),
                                        c(120,126), c(124,128,130))
Usatnts(testgen) <- c(2,2,2)

# look up which samples*loci have missing genotypes
find.missing.gen(testgen)

# get the three-dimensional distance array and the mean of the array
gendist <- meandistance.matrix(testgen, distmetric=Bruvo.distance,
                                maxl=4, all.distances=TRUE)
# look at the distances for loc1, where there is missing data and long genotypes
gendist[[1]]["loc1",,]

# look up all missing distances in the array
find.na.dist(gendist[[1]])

# look up just the missing distances that don't result from missing genotypes
find.na.dist.not.missing(testgen, gendist[[1]])

# Copy the array to edit the new copy
newDistArray <- gendist[[1]]
# calculate the distances that were NA from genotype lengths exceeding maxl
# (in reality, if this were too computationally intensive you might estimate
# it manually instead)
subDist <- Bruvo.distance(c(100,106,108,110,114), c(102,104,106,110,112))
subDist
# insert this distance into the correct positions
newDistArray["loc1","ind3","ind4"] <- subDist
newDistArray["loc1","ind4","ind3"] <- subDist
# calculate the new mean distance matrix
newMeanMatrix <- meandist.from.array(newDistArray)
# look at the difference between this matrix and the original.
newMeanMatrix
gendist[[2]]
}
\keyword{ NA }
\keyword{ array }
\keyword{ arith }
