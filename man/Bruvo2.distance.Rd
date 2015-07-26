\name{Bruvo2.distance}
\alias{Bruvo2.distance}
\title{
Distance Measure of Bruvo et al. under Genome Loss and Addition
}
\description{
This is an inter-individual distance measure similar to
\code{\link{Bruvo.distance}}, except that where genotypes have different
numbers of alleles, virtual alleles are equal to those from the longer
and/or shorter genotype, rather than being equal to infinity.
}
\usage{
Bruvo2.distance(genotype1, genotype2, maxl = 7, usatnt = 2,
                missing = -9, add = TRUE, loss = TRUE)
}
\arguments{
  \item{genotype1}{
A numeric vector representing the genotype of one individual at one
locus.  This type of vector is produced by the \code{\link{Genotype}}
method for \code{"genambig"} objects.
}
  \item{genotype2}{
The second genotype for the distance calculation, in the same format as
\code{genotype1}.
}
  \item{maxl}{
The maximum number of alleles that either genotype can have.  If it is
exceeded, \code{NA} is returned.  This argument exists to prevent
computations that would take in excess of an hour; see
\code{\link{Bruvo.distance}}.
}
  \item{usatnt}{
The microsatellite repeat type for the locus.  \code{2} for dinucleotide
repeats, \code{3} for dinucleotide repeats, \code{1} if the alleles are
already coded as repeat numbers, etc.  See \code{\link{Usatnts}}.
}
  \item{missing}{
The symbol that indicates missing data for a given sample and locus.
See \code{\link{Missing}}.
}
  \item{add}{
\code{TRUE} if the model of genome addition is being used, and
\code{FALSE} if it is not.  If this model is used, the shorter genotype
will have virtual alleles added from the same genotype.
}
  \item{loss}{
\code{TRUE} if the model of genome loss is being used, and \code{FALSE}
if it is not.  If this model is used, the shorter genotype will have
virtual alleles added from the longer genotype.
}
}
\details{
Bruvo \emph{et al.} (2004) describe multiple methods for calculating
genetic distances between individuals of different ploidy.  (See
\dQuote{Special cases} starting on page 2102 of the paper.)  The
original \code{Bruvo.distance} function in \pkg{polysat} uses the method
described for systems with complex changes in ploidy, adding virtual
alleles equal to infinity to the shorter genotype to make it the same
length as the longer genotype.  This method, however, can exaggerate
distances between individuals of different ploidy, particularly when
used with \code{\link{meandistance.matrix2}} as opposed to
\code{\link{meandistance.matrix}}.

\code{Bruvo2.distance} calculates distances between individuals under
the models of genome addition and genome loss.  If \code{add = TRUE} and
\code{loss = TRUE}, the distance produced is equal to that of equation 6
in the paper.

If \code{add = TRUE} and \code{loss = FALSE}, the distance calculated is
that under genome addition only.  Likewise if \code{add = FALSE} and
\code{loss = TRUE} the distance is calculated under genome loss only.
The latter distance should be greater than the former.  If both were averaged
together, they would give the identical result to that produced when
\code{add = TRUE} and \code{loss = TRUE}.  All three distances will be
less than that produced by \code{Bruvo.distance}.

If both genotypes have the same number of alleles, they are passed to
\code{Bruvo.distance} for the calculation.  This also happens if
\code{add = FALSE} and \code{loss = FALSE}.  Otherwise, if the genotypes have
different numbers of alleles, all possible genotypes with virtual
alleles are enumerated and passed to \code{Bruvo.distance} one by one, and
the results averaged.

The number of different genotypes simulated under genome loss or genome
addition is \eqn{l^d}, where \eqn{l} is the length of the genotype from
which virtual alleles are being taken, and \eqn{d} is the difference in
length between the longer and shorter genotype.  For example, under
genome addition for a diploid individual with alleles 1 and 2 being
compared to a tetraploid individual, the genotypes 1211, 1212, 1221, and
1222 will each be used once to represent the diploid individual.
}
\value{
A decimal between 0 and 1, with 0 indicating complete identity of two
genotypes, and 1 indicating maximum dissimilarity.  \code{NA} is
returned if one or both genotypes are missing or if \code{maxl} is exceeded.
}
\references{
  Bruvo, R., Michiels, N. K., D'Sousa, T. G., and Schulenberg, H. (2004)
  A simple method for calculation of microsatellite genotypes
  irrespective of ploidy level. \emph{Molecular Ecology} \bold{13}, 2101--2106.
}
\author{
Lindsay V. Clark
}
\note{
Figure 1B and 1C of Bruvo \emph{et al.} (2004) illustrate an example of
this distance measure.  To perform the identical calculation to that
listed directly under equation 6, you would type:

\code{Bruvo2.distance(c(20,23,24), c(20,24,26,43), usatnt=1)}

However, you will notice that the result, 0.401, is slightly different from that
given in the paper.  This is due to an error in the paper.  For the
distance under genome loss when the virtual allele is 26, the result
should be 1 instead of 1.75.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{Lynch.distance}}, \code{\link{Bruvo.distance}},
  \code{\link{meandistance.matrix2}}
}
\examples{
Bruvo2.distance(c(102,104), c(104,104,106,110))
Bruvo2.distance(c(102,104), c(104,104,106,110), add = FALSE)
Bruvo2.distance(c(102,104), c(104,104,106,110), loss = FALSE)
}
\keyword{ arith }

