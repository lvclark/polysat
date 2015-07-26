\name{simpleFreq}
\alias{simpleFreq}
\title{
Simple Allele Frequency Estimator
}
\description{
Given genetic data, allele frequencies by population are calculated.
This estimation method assumes polysomic inheritance.  For genotypes
with allele copy number ambiguity, all alleles are assumed to have an
equal chance of being present in multiple copies.  This function is best
used to generate initial values for more complex allele frequency
estimation methods.
}
\usage{
simpleFreq(object, samples = Samples(object), loci = Loci(object))
}
\arguments{
  \item{object}{
A \code{genbinary} or \code{genambig} object containing genotype data.
No \code{NA} values are allowed for \code{PopInfo(object)[samples]} or
\code{Ploidies(object,samples,loci)}. (Population identity and ploidy are
needed for allele frequency calculation.)
}
  \item{samples}{
An optional character vector of samples to include in the calculation.
}
  \item{loci}{
An optional character vector of loci to include in the calculation.
}
}
\details{
If \code{object} is of class \code{genambig}, it is converted to a
\code{genbinary} object before allele frequency calculations take
place.  Everything else being equal, the function will work more quickly
if it is supplied with a \code{genbinary} object.

  For
each sample*locus, a conversion factor is generated that is the ploidy
of the sample (and/or locus) as specified in \code{Ploidies(object)}
divided by the number of
alleles that the sample has at that locus.  Each allele is then
considered to be present in as many copies as the conversion factor
(note that this is not necessarily an integer).  The number of copies of
an allele is totaled for the population and is divided by the total
number of genomes in the population (minus missing data at the locus)
in order to calculate allele frequency.

A major assumption of this calculation method is that each allele in a
partially heterozygous genotype has an equal chance of being present in
more than one copy.  This is almost never true, because common alleles
in a population are more likely to be partially homozygous in an
individual.  The result is that the frequency of common alleles is
underestimated and the frequency of rare alleles is overestimated.  Also
note that the level of inbreeding in the population has an effect on the
relationship between genotype frequencies and allele frequencies, but is
not taken into account in this calculation.
}
\value{
Data frame, where each population is in one row.  If each sample in
\code{object} has only one ploidy, the first column of the data frame is
called \code{Genomes} and contains the number of genomes in each
population.  Otherwise, there is a
\code{Genomes} column for each locus.  Each remaining column contains
frequencies for one allele.
Columns are named by locus and allele, separated by a period.  Row names
are taken from \code{PopNames(object)}.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\linkS4class{genbinary}}, \code{\linkS4class{genambig}}
}
\examples{
# create a data set for this example
mygen <- new("genambig", samples = paste("ind", 1:6, sep=""),
             loci = c("loc1", "loc2"))
mygen <- reformatPloidies(mygen, output="sample")
Genotypes(mygen, loci="loc1") <- list(c(206),c(208,210),c(204,206,210),
    c(196,198,202,208),c(196,200),c(198,200,202,204))
Genotypes(mygen, loci="loc2") <- list(c(130,134),c(138,140),c(130,136,140),
    c(138),c(136,140),c(130,132,136))
PopInfo(mygen) <- c(1,1,1,2,2,2)
Ploidies(mygen) <- c(2,2,4,4,2,4)

# calculate allele frequencies
myfreq <- simpleFreq(mygen)

# look at the results
myfreq

# an example where ploidy is indexed by locus instead
mygen2 <- new("genambig", samples = paste("ind", 1:6, sep=""),
             loci = c("loc1", "loc2"))
mygen2 <- reformatPloidies(mygen2, output="locus")
PopInfo(mygen2) <- 1
Ploidies(mygen2) <- c(2,4)
Genotypes(mygen2, loci="loc1") <- list(c(198), c(200,204), c(200),
                                       c(198,202), c(200), c(202,204))
Genotypes(mygen2, loci="loc2") <- list(c(140,144,146), c(138,144),
                                       c(136,138,144,148), c(140),
                                       c(140,142,146,150),
                                       c(142,148,150))
myfreq2 <- simpleFreq(mygen2)
myfreq2
}

\keyword{ arith }

