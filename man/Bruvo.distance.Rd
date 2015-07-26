\name{Bruvo.distance}
\alias{Bruvo.distance}
\title{Genetic Distance Metric of Bruvo et al.}
\description{
  This function calculates the distance between two individuals at one
  microsatellite locus using a method based on that of Bruvo
  \emph{et al.} (2004).
}
\usage{Bruvo.distance(genotype1, genotype2, maxl=9, usatnt=2, missing=-9)}
\arguments{
  \item{genotype1}{A vector of alleles for one individual at one
     locus.  Allele length is in nucleotides or repeat count.  Each
     unique allele corresponds to one element in the vector, and the
     vector is no longer than it needs to be to contain all unique
     alleles for this individual at this locus.}
  \item{genotype2}{A vector of alleles for another individual at the
      same locus.}
  \item{maxl}{If both individuals have more than this number of
    alleles at this locus, \code{NA} is returned instead of a
    numerical distance.}
  \item{usatnt}{Length of the repeat at this locus.  For example
    \code{usatnt=2} for dinucleotide repeats, and \code{usatnt=3}
    for trinucleotide repeats.  If the alleles in \code{genotype1}
    and \code{genotype2} are expressed in repeat count instead of
    nucleotides, set \code{usatnt=1}.}
  \item{missing}{A numerical value that, when in the first allele
    position, indicates missing data. \code{NA} is returned if this
    value is found in either genotype.}
}
\details{
  Since allele copy number is frequently unknown in polyploid
  microsatellite data, Bruvo \emph{et al.} developed a measure of genetic
  distance similar to band-sharing indices used with dominant data,
  but taking into account mutational distances between alleles.  A
  matrix is created containing all differences in repeat count between
  the alleles of two individuals at one locus.  These differences are
  then geometrically transformed to reflect the probabilities of
  mutation from one allele to another.  The matrix is then searched to
  find the minimum sum if each allele from one individual is paired to
  one allele from the other individual.  This sum is divided by the
  number of alleles per individual.

  If one genotype has more alleles than the other, \sQuote{virtual alleles} must
  be created so that both genotypes are the same length.  There are
  three options for the value of these virtual alleles, but
  \code{Bruvo.distance} only implements the simplest one, assuming that it is
  not known whether differences in ploidy arose from genome addition or
  genome loss.  Virtual alleles are set to infinity, such that the
  geometric distance between any allele and a virtual allele is 1.

  In the original publication by Bruvo \emph{et al.} (2004), ambiguous
  genotypes were dealt with by calculating the distance for all possible
  unambiguous genotype combinations and averaging across all of them
  equally.  When \code{Bruvo.distance} is called from
  \code{meandistance.matrix}, ploidy is unknown and all genotypes are
  simply treated as if they had one copy of each allele.  When
  \code{Bruvo.distance} is called from \code{meandistance.matrix2}, the
  analysis is truer to the original, in that ploidy is known and all
  possible unambiguous genotype combinations are considered.  However,
  instead of all possible unambiguous genotypes being weighted equally,
  in \code{meandistance.matrix2} they are weighted based on allele
  frequencies and selfing rate, since some unambiguous genotypes are
  more likely than others.
}
\value{
  A number ranging from 0 to 1, with 0 indicating identical
    genotypes, and 1 being a theoretical maximum distance if all alleles
    from \code{genotype1} differed by an infinite number of repeats from all
    alleles in \code{genotype2}.  \code{NA} is returned if both genotypes have
    more than \code{maxl} alleles or if either genotype has the symbol for
    missing data as its first allele.
}
\references{
  Bruvo, R., Michiels, N. K., D'Sousa, T. G., and Schulenberg, H. (2004)
  A simple method for calculation of microsatellite genotypes
  irrespective of ploidy level. \emph{Molecular Ecology} \bold{13}, 2101-2106.
}
\note{
  The processing time is a function of the factorial of the number of
  alleles, since each possible combination of allele pairs must be
  evaluated.  For genotypes with a sufficiently large number of alleles,
  it may be more efficient to estimate distances manually by creating
  the matrix in Excel and visually picking out the shortest distances
  between alleles.  This is the purpose of the \code{maxl} argument.  On my
  personal computer, if both genotypes had more than nine alleles, the
  calculation could take an hour or more, and so this is the default
  limit.  In this case, \code{Bruvo.distance} returns \code{NA}.
}
\seealso{
  \code{\link{meandistance.matrix}}, \code{\link{Lynch.distance}},
  \code{\link{Bruvo2.distance}}
}
\examples{
  Bruvo.distance(c(202,206,210,220),c(204,206,216,222))
  Bruvo.distance(c(202,206,210,220),c(204,206,216,222),usatnt=4)
  Bruvo.distance(c(202,206,210,220),c(204,206,222))
  Bruvo.distance(c(202,206,210,220),c(204,206,216,222),maxl=3)
  Bruvo.distance(c(202,206,210,220),c(-9))
}
\author{Lindsay V. Clark}
\keyword{arith}
