\name{genotypeDiversity}
\alias{genotypeDiversity}
\alias{Shannon}
\alias{Simpson}
\alias{Simpson.var}
\title{
Genotype Diversity Statistics
}
\description{
\code{genotypeDiversity} calculates diversity statistics based on
genotype frequencies, using a distance matrix to assign individuals to
genotypes.  The \code{Shannon} and \code{Simpson} functions are also
available to calculate these statistics directly from a vector
of frequencies.
}
\usage{
genotypeDiversity(genobject, samples = Samples(genobject),
                  loci = Loci(genobject),
                  d = meandistance.matrix(genobject, samples, loci,
                                          all.distances = TRUE,
                                          distmetric = Lynch.distance),
                  threshold = 0, index = Shannon, ...)

Shannon(p, base = exp(1))

Simpson(p)

Simpson.var(p)
}
\arguments{
  \item{genobject}{
An object of the class \code{"genambig"} (or more generally,
\code{"gendata"} if a value is supplied to \code{d}).  If there is more
than one population, the \code{PopInfo} slot should be filled in.
\code{genobject}
is the dataset to be analyzed, although the genotypes themselves will
not be used if \code{d} has already been calculated.  Missing genotypes,
however, will indicate individuals that should be skipped in the analysis.
}
  \item{samples}{
An optional character vector indicating a subset of samples to analyze.
}
  \item{loci}{
An optional character vector indicating a subset of loci to analyze.
}
  \item{d}{
A list such as that produced by \code{\link{meandistance.matrix}} or
\code{\link{meandistance.matrix2}} when \code{all.distances = TRUE}.  The
first item in the list is a three dimensional array, with the first
dimension indexed by locus and the second and third dimensions indexed
by sample.  These are genetic distances between samples, by locus.  The
second item in the list is the distance matrix averaged across loci.
This mean matrix will be used only if all loci are being analyzed.  If
\code{loci} is a subset of the loci found in \code{d}, the mean matrix
will be recalculated.
}
  \item{threshold}{
The maximum genetic distance between two samples that can be considered
to be the same genotype.
}
  \item{index}{
The diversity index to calculate.  This should be \code{Shannon},
\code{Simpson}, or a user-defined function that takes as its first
argument a vector of frequencies that sum to one.
}
  \item{\dots}{
Additional arguments to pass to \code{index}, for example the
\code{base} argument for \code{Shannon}.
}
\item{p}{A vector of counts.}
\item{base}{
The base of the logarithm for calculating the Shannon index.  This is
\code{exp(1)} for the natural log, or \code{2} for log base 2.
}
}
\details{
\code{genotypeDiversity} runs \code{\link{assignClones}} on distance
matrices for individual loci and then for all loci, for each seperate
population.  The results of \code{assignClones} are used to
calculate a vector of genotype frequencies, which is passed to
\code{index}.

\code{Shannon} calculates the Shannon index, which is:

\deqn{-\sum \frac{p_i}{N}\ln(\frac{p_i}{N})}{-\sum p_i/N ln(p_i/N)}

(or log base 2 or any other base, using the \code{base} argument) given
a vector \eqn{p} of genotype counts, where \eqn{N} is the sum of those counts.

\code{Simpson} calculates the Simpson index, which is:

\deqn{\sum \frac{p_{i}(p_{i} - 1)}{N(N - 1)}}{\sum p_i(p_i - 1)/(N(N -1))}

\code{Simpson.var} calculates the variance of the Simpson index:

\deqn{\frac{4N(N-1)(N-2)\sum p_{i}^3 + 2N(N-1)\sum p_{i}^2 -
    2N(N-1)(2N-3)(\sum p_{i}^2)^2}{[N(N-1)]^2}}

The variance of the Simpson index can be used to calculate a confidence
interval, for example the results of \code{Simpson} plus or minus twice
the square root of the results of \code{Simpson.var} would be the 95\%
confidence interval.
}
\value{
A matrix of diversity index results, with populations in rows and
loci in columns.  The final column is called \code{"overall"} and gives
the results when all loci are analyzed together.
}
\references{
Shannon, C. E. (1948) A mathematical theory of communication.
\emph{Bell System Technical Journal} \bold{27}:379--423 and 623--656.

Simpson, E. H. (1949) Measurement of diversity.  \emph{Nature}
\bold{163}:688.

Lowe, A., Harris, S. and Ashton, P. (2004) \emph{Ecological Genetics:
  Design, Analysis, and Application}.  Wiley-Blackwell.

Arnaud-Haond, S., Duarte, M., Alberto, F. and Serrao, E. A.  (2007)
Standardizing methods to address clonality in population studies.
\emph{Molecular Ecology} \bold{16}:5115--5139.

\url{http://darwin.phyloviz.net/ComparingPartitions/index.php?link=Tut4}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{assignClones}}, \code{\link{alleleDiversity}}
}
\examples{
# set up dataset
mydata <- new("genambig", samples=c("a","b","c"), loci=c("F","G"))
Genotypes(mydata, loci="F") <- list(c(115,118,124),c(115,118,124),
                                   c(121,124))
Genotypes(mydata, loci="G") <- list(c(162,170,174),c(170,172),
                                    c(166,180,182))
Usatnts(mydata) <- c(3,2)

# get genetic distances
mydist <- meandistance.matrix(mydata, all.distances=TRUE)

# calculate diversity under various conditions
genotypeDiversity(mydata, d=mydist)
genotypeDiversity(mydata, d=mydist, base=2)
genotypeDiversity(mydata, d=mydist, threshold=0.3)
genotypeDiversity(mydata, d=mydist, index=Simpson)
genotypeDiversity(mydata, d=mydist, index=Simpson.var)
}
\keyword{ arith }

