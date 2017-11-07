\name{alleleCorrelations}
\alias{alleleCorrelations}
\alias{testAlGroups}
\title{
Assign Alleles to Isoloci Based on Distribution of Genotypes
}
\description{
Where a single locus represents two or more independent isoloci (as in
an allopolyploid, or a diploidized autopolyploid), these two functions
can be used in sequence to assign alleles to isoloci.
\code{alleleCorrelations} uses K-means and UPGMA clustering of pairwise p-values
from Fisher's exact test to make initial groupings of alleles into
putative isoloci.  \code{testAlGroups} is then used to check those
groupings against individual genotypes, and adjust the assignments if necessary.
}
\usage{
alleleCorrelations(object, samples = Samples(object), locus = 1,
                   alpha = 0.05, n.subgen = 2, n.start = 50)

testAlGroups(object, fisherResults, SGploidy=2, samples=Samples(object),
             null.weight=0.5, tolerance=0.05, swap = TRUE,
             R = 100, rho = 0.95, T0 = 1, maxreps = 100)
}
\arguments{
  \item{object}{
A \code{"\linkS4class{genambig}"} or \code{"\linkS4class{genbinary}"}
object containing the data to analyze.
}
  \item{samples}{
An optional character or numeric vector indicating which samples to analyze.
}
  \item{locus}{
A single character string or integer indicating which locus to analyze.
}
  \item{alpha}{
The significance threshold, before multiple correction, for determining
whether two alleles are significantly correlated.
}
  \item{n.subgen}{
The number of subgenomes (number of isoloci) for this locus.  This would
be \code{2} for an allotetraploid or \code{3]} for an allohexaploid.  For an 
allo-octoploid, the value would be \code{2} if there were two tetraploid 
subgenomes, or \code{4} if there were four diploid subgenomes.
}
  \item{n.start}{
Integer, passed directly to the \code{nstart} argument of the R base
function \code{kmeans}.  Lowering this number will speed up computation
time, whereas increasing it will improve the probability of finding the
correct allele assignments.  The default value of 50 should work well in
most cases.
}
\item{fisherResults}{
  A list output from \code{alleleCorrelations}.
}
\item{SGploidy}{
  The ploidy of each subgenome (each isolocus).  This is \code{2} for an
    allotetraploid, an allohexaploid, or an allo-octoploid with four tetraploid 
    subgenomes, or \code{4} for an allo-octoploid with
    two tetraploid genomes.
  }
  \item{null.weight}{
    Numeric, indicating how genotypes with potential null alleles should
    be counted when looking for signs of homoplasy.  \code{null.weight}
    should be \code{0} if null
    alleles are expected to be common, and \code{1} if there are no null
    alleles in the dataset.  The default of \code{0.5} was chosen to reflect
    the fact that the presence of null alleles is generally unknown.
  }
  \item{tolerance}{
    The proportion of genotypes that are allowed to be in disagreement
    with the allele assignments.  This is the proportion of genotypes
    that are expected to have meiotic error or scoring error.
  }
  \item{swap}{
    Boolean indicating whether or not to use the allele swapping algorithm
    before checking for homoplasy.  TRUE will yield more accurate results
    in most cases, but FALSE may be preferable for loci with null or 
    homoplasious alleles at high frequency.
  }
  \item{R}{
    Simulated annealing parameter for the allele swapping algorithm.
    Indicates how many swaps to attempt in each rep (\emph{i.e.} how many
    swaps to attempt before changing the temperature).
  }
  \item{rho}{
    Simulated annealing parameter for the allele swapping algorithm.
    Factor by which to reduce the temperature at the end of each rep.
  }
  \item{T0}{
    Simulated annealing parameter for the allele swapping algorithm.
    Starting temperature.
  }
  \item{maxreps}{
    Simulated annealing parameter for the allele swapping algorithm.
    Maximum number of reps if convergence is not achieved.
  }
}
\details{
These functions implement a novel methodology, introduced in
\pkg{polysat} version 1.4 and updated in version 1.6, for cases where 
one pair of microsatellite
primers amplifies alleles at two or more independently-segregating loci
(referred to here as isoloci).  This is not typically the case with new
autopolyploids, in which all copies of a locus have equal chances of
pairing with each other at meiosis.  It is, however, frequently the case
with allopolyploids, in which there are two homeologous subgenomes that
do not pair (or infrequently pair) at meiosis, or ancient
autopolyploids, in which duplicated chromosomes have diverged to the
point of no longer pairing at meiosis.

Within the two functions there are four major steps:

\enumerate{
  \item \code{alleleCorrelations} checks to see if there are any alleles
  that are present in every genotype in the dataset.  Such invariable
  alleles are assumed to be fixed at one isolocus (which is not
  necessarily true, but may be corrected by
  \code{testAlGroups} in steps 4 and 5).
  If present, each invariable allele is assigned to its own isolocus.
  If there are more invariable alleles than isoloci, the function throws
  an error.  If only one isolocus remains, all remaining (variable) alleles are
  assigned to that isolocus.  If there are as many invariable alleles as
  isoloci, all remaining (variable) alleles are assigned to all isoloci
  (\emph{i.e.} they are considered homoplasious because they cannot be
  assigned).
  \item If, after step 1, two or more isoloci remain
  without alleles assigned to them, correlations between alleles are
  tested by \code{alleleCorrelations}.  The dataset is converted
  to \code{"genbinary"} if not
  already in that format, and a Fisher's exact test, with negative
  association (odds ratio being less than one) as the alternative
  hypothesis, is performed between
  each pair of columns (alleles) in the genotype matrix.  The p-value of
  this test between each pair of alleles is stored in a square matrix,
  and zeros are inserted into the diagonal of the matrix.  K-means
  clustering and UPGMA are then performed on the square matrix of
  p-values, and the
  clusters that are produced represent initial assignments of alleles
  to isoloci.  
  \item The output of \code{alleleCorrelations} is then passed to
  \code{testAlGroups}.  If the results of K-means clustering and UPGMA
  were not identical, \code{testAlGroups} checks both sets of
  assignments against all genotypes in the dataset.  For a genotype to
  be consistent with a set of assignments, it should have at least one
  allele and no more than \code{SGploidy} alleles belonging to each
  isolocus.  The set of assignments that is consistent with the greatest
  number of genotypes is chosen, or in the case of a tie, the set of
  assignments produced by K-means clustering.
  \item If \code{swap = TRUE} and the assignments chosen in the previous 
  step are inconsistent with some genotypes, \code{testAlGroups} attempts
  to swap the isoloci of single alleles, using a simulated annealing 
  (Bertsimas and Tsitsiklis 1993) algorithm to search for a new set of 
  assignments that is consistent with as many genotypes as possible.
  At each step, an allele is chosen at random to be moved to a different
  isolocus (which is also chosen at random if there are more than two
  isoloci).  If the new set of allele assignments is consistent with an equal or 
  greater number of genotypes than the previous set of assignments, the new
  set is retained.  If the new set is consistent with fewer genotypes than
  the old set, there is a small probability of retaining the new set, 
  dependent on how much worse the new set of assignments is and what the
  current \dQuote{temperature} of the algorithm is.  After \code{R} allele
  swapping attempts, the temperature is lowered, reducing the probability 
  of retaining a set of allele assignments that is worse than the previous set.
  A new rep of \code{R} swapping attempts then begins.
  If a set of allele assignments is found that is consistent with all genotypes,
  the algorithm stops immediately.  Otherwise it stops if no changes are made
  during an entire rep of \code{R} swap attempts, or if \code{maxreps} reps
  are performed.
  \item \code{testAlGroups} then checks through all genotypes to look
  for signs of homoplasy, meaning single alleles that should be assigned
  to more than one isolocus.  For each genotype, there should be no more
  than \code{SGploidy} alleles assigned to each isolocus.  Additionally,
  if there are no null alleles, each genotype should have at least one
  allele belonging to each isolocus.  Each time a genotype is
  encountered that does not meet these criteria, the a score is
  increased for all alleles that might be homoplasious.  (The second
  criterion is not checked if \code{null.weight = 0}.)  This score
  starts at zero and is increased by 1 if there are too many alleles per
  isolocus or by \code{null.weight} if an isolocus has no alleles.  Once
  all genotypes have been checked, the allele with the highest score is
  considered to be homoplasious and is added to the other isolocus.  (In
  a hexaploid or higher, which isolocus the allele is added to depends on the
  genotypes that were found to be inconsistent with the allele
  assignments, and which isolocus or isoloci the allele could have
  belonged to in order to fix the assignment.)  Allele scores are reset
  to zero and all alleles are then
  checked again with the new set of allele assignments.  The process is
  repeated until the proportion of genotypes that are inconsistent with
  the allele assignments is at or below \code{tolerance}.
  }
}
\value{
  Both functions return lists.  For \code{alleleCorrelations}:
  \item{locus }{The name of the locus that was analyzed.}
  \item{clustering.method }{The method that was ultimately used to
    produce \code{value$Kmeans.groups} and \code{value$UPGMA.groups}.
    Either \code{"K-means and UPGMA"} or \code{"fixed alleles"}.}
  \item{significant.neg }{Square matrix of logical values indicating
  whether there was significant negative correlation between each pair
  of alleles, after multiple testing correction by Holm-Bonferroni.}
\item{significant.pos }{Square matrix of logical values indicating
  whether there was significant positive correlation between each pair
  of alleles, after multiple testing correction by Holm-Bonferroni.}
\item{p.values.neg }{Square matrix of p-values from Fisher's exact test
  for negative correlation between each pair of alleles.}
\item{p.values.pos }{Square matrix of p-values from Fisher's exact test
  for positive correlation between each pair of alleles.}
\item{odds.ratio }{Square matrix of the odds ratio estimate from
  Fisher's exact test for each pair of alleles.}
\item{Kmeans.groups }{Matrix with \code{n.subgen} rows, and as many
  columns as there are alleles in the dataset.  \code{1} indicates that
  a given allele belongs to a given isolocus, and \code{0} indicates
  that it does not.  These are the groupings determined by K-means
  clustering.}
\item{UPGMA.groups }{Matrix in the same format as
  \code{value$Kmeans.groups}, showing groupings determined by UPGMA.}
\item{heatmap.dist }{Square matrix like \code{value$p.values.neg} but
  with zeros inserted on the diagonal.  This is the matrix that was used
  for K-means clustering and UPGMA.  This matrix can be passed to the
  \code{heatmap} function in R to visualize the clusters.}
\item{totss }{Total sums of squares output from K-means clustering.}
\item{betweenss }{Sums of squares between clusters output from K-means
  clustering.  \code{value$betweenss/value$totss} can be used as an
  indication of clustering quality.}
\item{gentable }{The table indicating presence/absence of each allele in
  each genotype.}

For \code{testAlGroups}:
\item{locus }{Name of the locus that was tested.}
\item{SGploidy }{The ploidy of each subgenome, taken from the
  \code{SGploidy} argument that was passed to \code{testAlGroups}.}
\item{assignments }{Matrix with as many rows as there are isoloci, and as many
  columns as there are alleles in the dataset.  \code{1} indicates that
  a given allele belongs to a given isolocus, and \code{0} indicates
  that it does not.}
\item{proportion.inconsistent.genotypes }{A number ranging from zero to one 
indicating the proportion of genotypes from the dataset that are inconsistent
with \code{assignments}.}
}
\references{
Clark, L. V. and Drauch Schreier, A. (2017) Resolving microsatellite genotype ambiguity in
populations of allopolyploid and diploidized autopolyploid organisms
using negative correlations between allelic variables. \emph{Molecular Ecology Resources}, 
\bold{17}, 1090--1103.  DOI: 10.1111/1755-0998.12639.

Bertsimas, D. and Tsitsiklis, J.(1993) Simulated annealing. \emph{Statistical Science}
\bold{8}, 10--15.

}
\author{
Lindsay V. Clark
}
\note{
\code{alleleCorrelations} will print a warning to the console or to the
  standard output stream if a significant positive correlation is found
  between any pair of alleles.  (This is not a \dQuote{warning} in the
  technical sense usually used in R, because it can occur by random
  chance and I did not want it to cause \pkg{polysat} to fail package
  checks.)  You can see which allele pair(s) caused this warning by
  looking at \code{value$significant.pos}.  If you receive this warning for
  many loci, consider that there may be population structure in your
  dataset, and that you might split the dataset into multiple
  populations to test seperately.  If it happens at just a few loci,
  check to make sure there are not scoring problems such as stutter
  peaks being miscalled as alleles.  If it only happens at one locus and
  you can't find any evidence of scoring problems, two alleles may have
  been positively correlated simply from random chance, and the warning
  can be ignored.
  
  \code{alleleCorrelations} can also produce an actual warning stating
  \dQuote{Quick-TRANSfer stage steps exceeded maximum}.  This warning
  is produced internally by \code{\link{kmeans}} and may occur if many
  genotypes are similar, as in mapping populations.  It can be safely
  ignored.
}


\seealso{
  \code{\link{recodeAllopoly}}, \code{\link{mergeAlleleAssignments}},
  \code{\link{catalanAlleles}}, \code{\link{processDatasetAllo}}
}
\examples{
# randomly generate example data for an allotetraploid
mydata <- simAllopoly(n.alleles=c(5,5), n.homoplasy=1)
viewGenotypes(mydata)

# test allele correlations
# n.start is lowered in this example to speed up computation time
myCorr <- alleleCorrelations(mydata, n.subgen=2, n.start=10)
myCorr$Kmeans.groups
myCorr$clustering.method
if(!is.null(myCorr$heatmap.dist)) heatmap(myCorr$heatmap.dist)

# check individual genotypes 
# (low maxreps used in order to speed processing time for this example)
myRes <- testAlGroups(mydata, myCorr, SGploidy=2, maxreps = 5)
myRes$assignments
myRes2 <- testAlGroups(mydata, myCorr, SGploidy=2, swap = FALSE)
myRes2$assignments
}
\keyword{ cluster }
\keyword{ logic }
