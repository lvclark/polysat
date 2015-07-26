\name{catalanAlleles}
\alias{catalanAlleles}
\title{
Sort Alleles into Isoloci
}
\description{
\code{catalanAlleles} uses genotypes present in a \code{"genambig"}
object to sort alleles from one locus into two or more isoloci in an
allopolyploid or diploidized autopolyploid species.  Alleles are
determined to belong to different
isoloci if they are both present in a fully homozygous genotype.  If
necessary, heterozygous genotypes are also examined to resolve remaining
alleles.
}
\usage{
catalanAlleles(object, samples = Samples(object), locus = 1,
               n.subgen = 2, SGploidy = 2, verbose = FALSE)
}
\arguments{
  \item{object}{
A \code{"genambig"} object containing the dataset to analyze.  All
individuals should be the same ploidy, although the function does not
access the \code{Ploidies} slot.  Missing data are allowed.  For the
locus to be examined, no genotype should have fewer than \code{n.subgen}
alleles or more than \code{n.subgen*SGploidy} alleles.
}
  \item{samples}{
Optional argument indicating samples to be analyzed.  Can be integer or
character, as with other \pkg{polysat} functions.
}
  \item{locus}{
An integer or character string indicating which locus to analyze.
Cannot be a vector greater than length 1.  (The function will only
analyze one locus at a time.)
}
  \item{n.subgen}{
    The number of isoloci (number of subgenomes).  For example, \code{2}
    for an allotetraploid,
and \code{3} for an allohexaploid (three diploid genomes).
}
  \item{SGploidy}{
The ploidy of each genome.  Only one value is allowed; all genomes
must be the same ploidy.  \code{2} indicates that each subgenome is
diploid (as in an allotetraploid, or an allohexaploid with three diploid
genomes).
}
  \item{verbose}{
Boolean.  Indicates whether results, and if applicable, problematic
genotypes, should be printed to the console.
}
}
\details{
\code{catalanAlleles} implements and extends an approach used by Catalan
\emph{et al.} (2006) that sorts alleles from a duplicated microsatellite
locus into two or more isoloci (homeologous loci on different
subgenomes).  First, fully homozygous genotypes are identified and used
in analysis.  If a genotype has as many alleles as there are subgenomes
(for example, a genotype with two alleles in an allotetraploid species),
it is assumed to be fully homozygous and the alleles are assumed to
belong to different subgenomes.  If some alleles remain unassigned after
examination of all fully homozygous genotypes, heterozygous genotypes
are also examined to attempt to assign those remaining alleles.

For example, in an allotetraploid, if a genotype contains one unassigned
allele, and all other
alleles in the genotype are known to belong to one isolocus, the
unassigned allele can be assigned to the other isolocus.  Or, if two
alleles in a genotype belong to one isolocus, one allele belongs to the
other isolocus, and one allele is unassigned, the unassigned allele can
be assigned to the latter isolocus.  The function follows such logic
(which can be extended to higher ploidies) until all alleles can be
assigned, or returns a text string saying that
the allele assignments were unresolvable.

It is important to note that this method assumes no null alleles and no
homoplasy across isoloci.  If the function encounters evidence of
either it will not return allele assignments.  Null alleles and
homoplasy are real possibilities in any dataset, which means that this
method simply will not work for some microsatellite loci.

(Null alleles are those that do not produce a PCR amplicon, usually
because of a mutation in the primer binding site.  Alleles that exhibit
homoplasy are those that produce amplicons of the same size, despite not
being identical by descent.  Specifically, homoplasy between alleles from
different isoloci will interfere with the Catalan method of allele
assignment.)
}
\value{
A list containing the following items:
  \item{locus }{A character string giving the name of the locus.}
  \item{SGploidy }{A number giving the ploidy of each subgenome.
    Identical to the \code{SGploidy} argument.}
  \item{assignments }{If assignments cannot be made, a character string
    describing the problem.  Otherwise, a matrix with \code{n.subgen}
    rows and a labeled column for each allele, with a \code{1} if the
    allele belongs to that subgenome and a \code{0} if it does not.}
}
\references{
Catalan, P., Segarra-Moragues, J. G., Palop-Esteban, M., Moreno, C. and
Gonzalez-Candelas, F. (2006) A Bayesian approach for discriminating
among alternative inheritance hypotheses in plant polyploids: the
allotetraploid origin of genus Borderea (Dioscoreaceae).
\emph{Genetics} \bold{172}, 1939--1953.
}
\author{
Lindsay V. Clark
}
\note{
  Aside from homoplasy and null alleles, stochastic effects may prevent
  the minimum combination of genotypes needed to resolve all alleles
  from being present in the dataset.  For a typical allotetraploid
  dataset, 50 to 100 samples will be needed, whereas an
  allohexaploid dataset may require over 100 samples.  In simulations,
  allo-octoploid datasets with two tetraploid genomes were unresolvable
  even with 10,000 samples due to the low probability of finding full
  homozygotes.  Additionally, loci are less likely to be resolvable if
  they have many alleles or if one isolocus is monomorphic.
  
  Although determination of allele copy number by is not needed (or
  expected) for \code{catalanAlleles} as it was in the originally
  published Catalan method, it is still very important that the
  genotypes be high quality.  Even a single scoring error can cause the
  method to fail, including allelic dropout, contamination between
  samples, stutter peaks miscalled as alleles, and PCR artifacts
  miscalled as alleles.  Poor quality loci (those that require some
  \dQuote{artistic} interpretation of gels or electropherograms) are
  unlikely to work with this method.  Individual genotypes that are of
  questionable quality should be discarded before running the function.
}


\seealso{
  \code{\link{alleleCorrelations}}, \code{\link{mergeAlleleAssignments}},
  \code{\link{recodeAllopoly}}, \code{\link{simAllopoly}}
}
\examples{
# make the default simulated allotetraploid dataset
mydata <- simAllopoly()

# resolve the alleles
myassign <- catalanAlleles(mydata)
}
\keyword{ misc }

