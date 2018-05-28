\name{meandistance.matrix}
\alias{meandistance.matrix}
\alias{meandistance.matrix2}
\title{Mean Pairwise Distance Matrix}
\description{
  Given a \code{genambig} object, \code{meandistance.matrix}
  produces a symmetrical matrix of
  pairwise distances between samples, averaged across all loci.  An
  array of all distances prior to averaging may also be produced.
}
\usage{
meandistance.matrix(object, samples = Samples(object),
                    loci = Loci(object), all.distances=FALSE,
                    distmetric = Bruvo.distance, progress = TRUE,
                    ...)
meandistance.matrix2(object, samples = Samples(object),
                     loci = Loci(object),
                     freq = simpleFreq(object, samples, loci), self = 0,
                     all.distances = FALSE, distmetric = Bruvo.distance,
                     progress = TRUE, ...)
}
\arguments{
  \item{object}{A \code{genambig} object containing the genotypes to be
    analyzed.  If \code{distmetric = Bruvo.distance}, the \code{Usatnts}
    slot should be filled in.  For \code{meandistance.matrix2},
    \code{Ploidies} and \code{PopInfo} are also required.}
  \item{samples}{A character vector of samples to be analyzed.  These
    should be all or a subset of the sample names used in \code{object}.}
  \item{loci}{A character vector of loci to be analyzed.  These should
    be all or a subset of the loci names used in \code{object}.}
  \item{freq}{A data frame of allele frequencies such as that produced
    by \code{\link{simpleFreq}} or \code{\link{deSilvaFreq}}.}
  \item{self}{A number ranging from 0 to 1, indicating the rate of selfing.}
  \item{all.distances}{If \code{FALSE}, only the mean distance matrix will be
    returned.  If \code{TRUE}, a list will be returned containing an array of
    all distances by locus and sample as well as the mean distance matrix.}
  \item{distmetric}{The function to be used to calculate distances
    between genotypes.  \code{Bruvo.distance}, \code{Lynch.distance}, or
  a distance function written by the user.}
  \item{progress}{If \code{TRUE}, loci and samples will be printed to
    the console as distances are calculated, so that the user can
    monitor the progress of the computation.}
  \item{...}{Additional arguments (such as \code{maxl}, \code{add}, and
    \code{loss}) to pass to \code{distmetric}.}
}
\value{
  A symmetrical matrix containing pairwise distances between all
  samples, averaged across all loci.  Row and column names of the matrix
  will be the sample names provided in the \code{samples} argument.  If
  \code{all.distances=TRUE}, a list will be produced containing the above
  matrix as well as a three-dimensional array containing all distances
  by locus and sample.  The array is the first item in the list, and the
  mean matrix is the second.
}
\details{
  Each distance for the three-dimensional array is calculated only once,
  to save computation time.  Since the array (and resulting mean matrix)
  is symmetrical, the distance is written to two positions in the array
  at once.

  \code{meandistance.matrix} uses ambiguous genotypes exactly as they
  are, whereas \code{meandistance.matrix2} uses
  \code{\link{genotypeProbs}} to calculate all possible unambiguous
  genotypes and their probabilities under random mating or partial
  selfing.  The distance
  between each possible pair of unambiguous genotypes for the two
  samples is calculated with \code{distmetric} and weighted by the
  product of the probabilities of the two gentoypes.  As you might
  expect, \code{meandistance.matrix2} takes longer to process a given
  \code{"genambig"} object than \code{meandistance.matrix} does.
  Additionally, the distance between two identical ambiguous genotypes
  will be zero when calculated with \code{meandistance.matrix}, and
  greater than zero when calculated with \code{meandistance.matrix2},
  due to potential differences in copy number of the alleles.

  When \code{Bruvo.distance} is used,
  \code{meandistance.matrix2} exaggerates distances between individuals
  of different ploidy as compared to \code{meandistance.matrix}.  The
  use of \code{Bruvo2.distance} with \code{meandistance.matrix2}
  allows individuals with different ploidies to have similar
  inter-individual distances to those between individuals of the same ploidy.
  In general, it will be desirable to use \code{Bruvo.distance} with
  \code{meandistance.matrix} for complex datasets with high ploidy
  levels, or \code{Bruvo.distance2} with \code{meandistance.matrix2} for
  hexaploid or lower datasets (based on how long it takes my personal
  computer to perform these calculations) where changes in ploidy are
  due to genome doubling or genome loss.  If all individuals have the
  same ploidy, \code{Bruvo.distance} and \code{Bruvo2.distance} will
  give identical results regardless of whether
  \code{meandistance.matrix} or \code{meandistance.matrix2} is used.

  \code{meandistance.matrix2} does not allow a genotype to have more
  alleles than the ploidy of the individual (as listed in the
  \code{Ploidies} slot).  Additionally, if \code{self} is greater than
  zero, each population may only have one ploidy at each locus.
}

\seealso{
  \code{\link{Bruvo.distance}}, \code{\link{Bruvo2.distance}},
  \code{\link{Lynch.distance}}, \code{\link{meandist.from.array}},
  \code{\link{GENLIST}}
}
\examples{
# create a list of genotype data
mygendata <- new("genambig", samples = c("ind1","ind2","ind3","ind4"),
                 loci = c("locus1","locus2","locus3","locus4"))
Genotypes(mygendata) <-
  array(list(c(124,128,138),c(122,130,140,142),c(122,132,136),c(122,134,140),
             c(203,212,218),c(197,206,221),c(215),c(200,218),
             c(140,144,148,150),c(-9),c(146,150),c(152,154,158),
             c(233,236,280),c(-9),c(-9),c(-9)))
Usatnts(mygendata) <- c(2,3,2,1)

# make index vectors of data to use
myloci <- c("locus1","locus2","locus3")
mysamples <- c("ind1","ind2","ind4")

# calculate array and matrix
mymat <- meandistance.matrix(mygendata, mysamples, myloci,
                             all.distances=TRUE)
# view the results
mymat[[1]]["locus1",,]
mymat[[1]]["locus2",,]
mymat[[1]]["locus3",,]
mymat[[2]]

# add addtional info needed for meandistance.matrix2
mygendata <- reformatPloidies(mygendata, output="one")
Ploidies(mygendata) <- 4
PopInfo(mygendata) <- c(1,1,1,1)

# calculate distances taking allele freqs into account
mymat2 <- meandistance.matrix2(mygendata, mysamples, myloci)
mymat2
# now do the same under selfing
mymat3 <- meandistance.matrix2(mygendata, mysamples, myloci, self=0.3)
mymat3
}
\author{Lindsay V. Clark}
\keyword{array}
\keyword{arith}
