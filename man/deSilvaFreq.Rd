\name{deSilvaFreq}
\alias{deSilvaFreq}
\title{
Estimate Allele Frequencies with EM Algorithm
}
\description{
This function uses the method of De Silva \emph{et al.} (2005) to
estimate allele frequencies under polysomic inheritance with a known
selfing rate.
}
\usage{
deSilvaFreq(object, self, samples = Samples(object),
            loci = Loci(object), initNull = 0.15,
            initFreq = simpleFreq(object[samples, loci]),
            tol = 1e-08, maxiter = 1e4)
}
\arguments{
  \item{object}{
    A \code{"\linkS4class{genambig}"} or \code{"\linkS4class{genbinary}"}
    object containing the dataset of
interest.  All ploidies for \code{samples} and \code{loci} should be the
same, and this should be an
even number.  \code{PopInfo} must also be filled in for \code{samples}.
}
  \item{self}{
A number between 1 and 0, indicating the rate of selfing.
}
  \item{samples}{
An optional character vector indicating a subset of samples to use in
the calculation.
}
  \item{loci}{
An optional character vector indicating a subset of loci for which to
calculate allele frequencies.
}
  \item{initNull}{
A single value or numeric vector indicating initial frequencies to use
for the null allele at each locus.
}
  \item{initFreq}{
A data frame containing allele frequencies (for non-null loci) to use
for initialization.  This needs to be in the same format as the output
of \code{\link{simpleFreq}} with a single \dQuote{Genomes} column
(similarly to the format of the
output of \code{deSilvaFreq}).  By default, the function will do a quick
estimation of allele frequencies using \code{simpleFreq} and then
initialize the EM algorithm at these frequencies.
}
  \item{tol}{
The tolerance level for determining when the results have converged.
Where \code{p2} and \code{p1} are the current and previous vectors of allele
frequencies, respectively, the EM algorithm stops if
\code{sum(abs(p2-p1)/(p2+p1)) <= tol}.
}
\item{maxiter}{
The maximum number of iterations that will be performed for each locus and 
population.
}
}
\details{
Most of the SAS code from the supplementary material of De Silva
\emph{et al.}
(2005) is translated directly into the R code for this function.  The
SIMSAMPLE (or CreateRandomSample in the SAS code) function is omitted
so that the actual allelic phenotypes from the dataset can be used
instead of simulated phenotypes.  \code{deSilvaFreq}
loops through each locus and population, and in each loop tallies the
number of alleles and sets up matrices using GENLIST, PHENLIST, RANMUL, SELFMAT,
and CONVMAT as described in the paper.
Frequencies of each allelic phenotype are then tallied
across all samples in that population with non-missing data at the
locus.  Initial allele
frequencies for that population and locus are then extraced from
\code{initFreq} and adjusted according to \code{initNull}.  The EM
iteration then begins for that population and locus, as described in the
paper (EXPECTATION, GPROBS, and MAXIMISATION).

Each repetition of the EM algorithm includes an expectation and
maximization step.  The expectation step uses allele frequencies and the
selfing rate to
calculate expected genotype frequencies, then uses observed phenotype
frequencies and expected genotype frequencies to estimate genotype
frequencies for the population.  The maximization step uses the
estimated genotype frequencies to calculate a new set of allele
frequencies.  The process is repeated until allele frequencies
converge.

In addition to returning a data frame of allele frequencies,
\code{deSilvaFreq} also prints to the console the number of EM
repetitions used for each population and locus.  When each locus and
each population is begun, a message is printed to the console so that
the user can monitor the progress of the computation.

}
\value{
A data frame containing the estimated allele frequencies.  The row names
are population names from \code{PopNames(object)}.  The first column
shows how many genomes each population has.  All other columns represent
alleles (including one null allele per locus).  These column names are
the locus name and allele name separated by a period.
}
\references{
  De Silva, H. N., Hall, A. J., Rikkerink, E., and Fraser, L. G. (2005)
  Estimation of allele frequencies in polyploids under
certain patterns of inheritance.  \emph{Heredity} \bold{95}, 327--334
}
\author{
Lindsay V. Clark
}
\note{
It is possible to exceed memory limits for R if a locus has too many
alleles in a population (e.g. 15 alleles in a tetraploid if the memory
limit is 1535 Mb, see \code{memory.limit}).

  De Silva \emph{et al.} mention that their estimation method could be
extended to the case of disomic inheritence.  A method for disomic
inheritence is not implemented here, as it would require knowledge of
which alleles belong to which isoloci.

De Silva \emph{et al.} also suggest a means of
estimating the selfing rate with a least-squares method.  Using the
notation in the source code, this would be:

\code{lsq <- smatt \%*\% EP - rvec}

\code{self <- as.vector((t(EP - rvec) \%*\% lsq)/(t(lsq) \%*\% lsq))}

However, in my experimentation with this calculation, it sometimes
	       yields selfing rates greater than one.  For this reason,
	       it is not implemented here.
	       }


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{simpleFreq}}, \code{\link{write.freq.SPAGeDi}},
  \code{\link{.genlist}}
}
\examples{
\dontrun{
## An example with a long run time due to the number of alleles

# create a dataset for this example
mygen <- new("genambig", samples=c(paste("A", 1:100, sep=""),
                                   paste("B", 1:100, sep="")),
             loci=c("loc1", "loc2"))
PopNames(mygen) <- c("PopA", "PopB")
PopInfo(mygen) <- c(rep(1, 100), rep(2, 100))
mygen <- reformatPloidies(mygen, output="one")
Ploidies(mygen) <- 4
Usatnts(mygen) <- c(2, 2)
Description(mygen) <- "An example for allele frequency calculation."

# create some genotypes at random for this example
for(s in Samples(mygen)){
    Genotype(mygen, s, "loc1") <- sample(seq(120, 140, by=2),
                                         sample(1:4, 1))
}
for(s in Samples(mygen)){
    Genotype(mygen, s, "loc2") <- sample(seq(130, 156, by=2),
                                         sample(1:4, 1))
}
# make one genotype missing
Genotype(mygen, "B4", "loc2") <- Missing(mygen)

# view the dataset
summary(mygen)
viewGenotypes(mygen)

# calculate the allele frequencies if the rate of selfing is 0.2
myfrequencies <- deSilvaFreq(mygen, self=0.2)

# view the results
myfrequencies
}

## An example with a shorter run time, for checking that the funciton
## is working.  Genotype simulation is also a bit more realistic here.

# Create a dataset for the example.
mygen <- new("genambig", samples=paste("A", 1:100, sep=""), loci="loc1")
PopNames(mygen) <- "PopA"
PopInfo(mygen) <- rep(1, 100)
mygen <- reformatPloidies(mygen, output="one")
Ploidies(mygen) <- 4
Usatnts(mygen) <- 2
for(s in Samples(mygen)){
    alleles <- unique(sample(c(122,124,126,0), 4, replace=TRUE,
                             prob = c(0.3, 0.2, 0.4, 0.1)))
    Genotype(mygen, s, "loc1") <- alleles[alleles != 0]
    if(length(Genotype(mygen, s, "loc1"))==0)
        Genotype(mygen, s, "loc1") <- Missing(mygen)
}

# We have created a random mating populations with four alleles
# including one null.  The allele frequencies are given in the
# 'prob' argument.

# Estimate allele frequencies
myfreq <- deSilvaFreq(mygen, self=0.01)
myfreq
}

\keyword{ array }
\keyword{ iteration }
