\name{assignClones}
\alias{assignClones}
\title{
Group Individuals Based on a Distance Threshold
}
\description{
\code{assignClones} uses a distance matrix such as that produced by
\code{\link{meandistance.matrix}} or \code{\link{meandistance.matrix2}}
to place individuals into groups representing asexually-related ramets,
or any other grouping based on a distance threshold.
}
\usage{
assignClones(d, samples = dimnames(d)[[1]], threshold = 0)
}

\arguments{
  \item{d}{
A symmetrical, square matrix containing genetic distances between
  individuals.  Both dimensions should be named according to the names
  of the individuals (samples).  A matrix produced by
  \code{meandistance.matrix} or \code{meandistance.matrix2} when
  \code{all.distances = FALSE}, or the matrix that is the second item in
  the list produced if \code{all.distances = TRUE}, will be in the right
  format. \code{\link{meandist.from.array}} will also produce a matrix
  in the correct format.
}
  \item{samples}{
A character vector containing the names of samples to analyze.  This
  should be all or a subset of the names of \code{d}.
}
  \item{threshold}{
A number indicating the maximum distance between two individuals that
  will be placed into the same group.
}
}
\details{
  This function groups individuals very similarly to the software GenoType
  (Meirmans and van Tienderen, 2004).  If a distance matrix from
  \pkg{polysat} is exported to GenoType, the results will be the same as
  those from \code{assignClones} assuming the same threshold is used.  Note
  that GenoType requires that
distances be integers rather than decimals, so you will have to
multiply the distances produced by \pkg{polysat} by a large number and
round them to the nearest integer if you wish to export them to GenoType.
When comparing the
results of \code{assignClones} and GenoType using my own data, the only
differences I have seen have been the result of rounding; a decimal that
was slightly above the threshold in when analyzed in R was rounded down
to the threshold when analyzed in GenoType.

Note that when using a distance threshold of zero (the default), it is
advisable to exclude all samples with missing data, in order to prevent
the merging of non-identical clones.  At higher thresholds, some missing
data are allowable, but samples that have missing data at many loci should
be excluded.

The \code{write.table} function can be used for exporting the results to
GenoDive.  See the R documentation for information on how to
make a tab-delimited file with no header.
}
\value{
A numeric vector, named by \code{samples}.  Each clone or group is given
  a number, and the number for each sample indicates the clone or group
  to which it belongs.
}
\references{
Meirmans, P. G. and Van Tienderen, P. H. (2004) GENOTYPE and GENODIVE:
  two programs for the analysis of genetic diversity of asexual
  organisms. \emph{Molecular Ecology Notes} \bold{4}, 792--794.
}
\author{
Lindsay V. Clark
}


\seealso{
\code{\link{genotypeDiversity}}
}
\examples{
# set up a simple matrix with three samples
test <- matrix(c(0,0,.5,0,0,.5,.5,.5,0), ncol=3, nrow=3)
abc <- c("a", "b", "c")
dimnames(test) <- list(abc,abc)

# assign clones with a threshold of zero or 0.5
assignClones(test)
assignClones(test, threshold=0.5)
}
\keyword{ arith }
