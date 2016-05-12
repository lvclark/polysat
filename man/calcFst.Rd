\name{calcFst}
\alias{calcFst}
\alias{calcPopDiff}
\title{
Estimate Population Differentiation Statistics
}
\description{
Given a data frame of allele frequencies and population sizes,
\code{calcPopDiff} calculates a matrix of pairwise \eqn{F_{ST}}{FST}, 
\eqn{G_{ST}}{GST}, or Jost's \eqn{D}{D} values.
\code{calcFst} is a wrapper for \code{calcPopDiff} to allow backwards
compatibility with previous versions of \pkg{polysat}.
}
\usage{
calcPopDiff(freqs, metric, pops = row.names(freqs), 
            loci = unique(as.matrix(
            as.data.frame(strsplit(names(freqs), split = ".", fixed = TRUE),
             stringsAsFactors = FALSE))[1, ]))

calcFst(freqs, pops = row.names(freqs), loci = unique(as.matrix(
        as.data.frame(strsplit(names(freqs), split = ".", fixed = TRUE),
        stringsAsFactors = FALSE))[1, ]))
}
\arguments{
  \item{freqs}{
A data frame of allele frequencies and population sizes such as that
produced by \code{simpleFreq} or \code{deSilvaFreq}.  Each population 
is in one row, and a
column called \code{Genomes} (or multiple columns containing the locus
names and \dQuote{Genomes} seperated by a period) contains the relative
size of each
population.  All other columns contain allele frequencies.  The names of
these columns are the locus name and allele name, separated by a period.
}
\item{metric}{
The population differentiation statistic to estimate.  Can be \dQuote{Fst},
\dQuote{Gst}, or \dQuote{Jost's D}.
}
  \item{pops}{
A character vector.  Populations to analyze, which should be
a subset of \code{row.names(freqs)}.
}
  \item{loci}{
A character vector indicating which loci to analyze.  These should be a
subset of the locus names as used in the column names of \code{freqs}.
}
}
\details{
For \code{metric = "Fst"} or \code{calcFst}:

\eqn{H_S}{HS} and \eqn{H_T}{HT} are estimate directly from allele frequencies
for each locus for each pair of populations, then averaged across loci.  
Wright's \eqn{F_{ST}}{FST} is then calculated for each pair of populations as 
\eqn{\frac{H_T - H_S}{H_T}}{(HT-HS)/HT}.

\eqn{H}{H} values (expected heterozygosities for populations and combined
populations) are calculated as one minus the sum of all squared allele
frequencies at a locus.  To calculte \eqn{H_T}{HT}, allele frequencies between two
populations are averaged before the calculation.  To calculate \eqn{H_S}{HS}, \eqn{H}{H}
values are averaged after the calculation.  In both cases, the averages
are weighted by the relative sizes of the two populations (as indicated
by \code{freqs$Genomes}).

For \code{metric = "Gst"}:

This metric is similar to \eqn{F_{ST}}{FST}, but mean allele frequencies and 
mean expected heterozygosities are not weighted by population size.  Additionally,
unbiased estimators of \eqn{H_S}{HS} and \eqn{H_T}{HT} are used according to Nei
and Chesser (1983; equations 15 and 16, reproduced also in Jost (2008)).
Instead of using twice the harmonic mean of the number of individuals in the two
subpopulations as described by Nei and Chesser (1983), the harmonic mean of the
number of allele copies in the two subpopulations (taken from \code{freq$Genomes})
is used, in order to allow for polyploidy.
\eqn{G_{ST}}{GST} is estimated for each locus and then averaged across loci.

For \code{metric = "Jost's D"}:

The unbiased estimators of \eqn{H_S}{HS} and \eqn{H_T}{HT} and calculated as 
with \eqn{G_{ST}}{GST}, without weighting by population size.  They are then used 
to estimate \eqn{D}{D} at each locus according to Jost (2008; equation 12):

\deqn{2 * \frac{H_T - H_S}{1 - H_S}}{2 * (HT - HS)/(1 - HS)}

Valued of \eqn{D}{D} are then averaged across loci.

}
\value{
A square matrix containing \eqn{F_{ST}}{FST}, \eqn{G_{ST}}{GST}, or \eqn{D}{D} values.  
The rows and columns of the matrix are both named by population.
}
\references{
Nei, M. (1973) Analysis of gene diversity in subdivided populations.
\emph{Proceedings of the National Academy of Sciences of the United
  States of America} \bold{70}, 3321--3323.
  
Nei, M. and Chesser, R. (1983) Estimation of fixation indices and gene diversities.
\emph{Annals of Human Genetics} \bold{47}, 253--259.
  
Jost, L. (2008) \eqn{G_{ST}}{GST} and its relatives to not measure differentiation.
\emph{Molecular Ecology} \bold{17}, 4015--4026.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{simpleFreq}}, \code{\link{deSilvaFreq}}
}
\examples{
# create a data set (typically done by reading files)
mygenotypes <- new("genambig", samples = paste("ind", 1:6, sep=""),
                   loci = c("loc1", "loc2"))
Genotypes(mygenotypes, loci = "loc1") <- list(c(206), c(208,210),
                                              c(204,206,210),
    c(196,198,202,208), c(196,200), c(198,200,202,204))
Genotypes(mygenotypes, loci = "loc2") <- list(c(130,134), c(138,140),
                                              c(130,136,140),
    c(138), c(136,140), c(130,132,136))
PopInfo(mygenotypes) <- c(1,1,1,2,2,2)
mygenotypes <- reformatPloidies(mygenotypes, output="sample")
Ploidies(mygenotypes) <- c(2,2,4,4,2,4)

# calculate allele frequencies
myfreq <- simpleFreq(mygenotypes)

# calculate pairwise differentiation statistics
myfst <- calcPopDiff(myfreq, metric = "Fst")
mygst <- calcPopDiff(myfreq, metric = "Gst")
myD <- calcPopDiff(myfreq, metric = "Jost's D")

# examine the results
myfst
mygst
myD
}
\keyword{ arith }
\keyword{ array }
