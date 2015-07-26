\name{calcFst}
\alias{calcFst}
\title{
Calculate Wright's Pairwise FST
}
\description{
Given a data frame of allele frequencies and population sizes,
\code{calcFst} calculates a matrix of pairwise Fst values.
}
\usage{
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
\code{calcFst} works by calculating HS and HT for each locus for each pair
of populations, then averaging HS and HT across loci.  FST is then
calculated for each pair of populations as (HT-HS)/HT.

H values (expected heterozygosities for populations and combined
populations) are calculated as one minus the sum of all squared allele
frequencies at a locus.  To calculte HT, allele frequencies between two
populations are averaged before the calculation.  To calculate HS, H
values are averaged after the calculation.  In both cases, the averages
are weighted by the relative sizes of the two populations (as indicated
by \code{freqs$Genomes}).
}
\value{
A square matrix containing FST values.  The rows and columns of the
matrix are both named by population.
}
\references{
Nei, M. (1973) Analysis of gene diversity in subdivided populations.
\emph{Proceedings of the National Academy of Sciences of the United
  States of America} \bold{70}, 3321--3323.
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

# calculate pairwise FST
myfst <- calcFst(myfreq)

# examine the results
myfst
}
\keyword{ arith }
\keyword{ array }
