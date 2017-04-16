\name{PIC}
\alias{PIC}
\title{
Polymorphic Information Content
}
\description{
Given a set of allele frequencies, this function estimates the Polymorphic Information Content (PIC) 
for each locus, within and/or across populations.
}
\usage{
PIC(freqs, pops = row.names(freqs), loci = unique(as.matrix(as.data.frame(
    strsplit(names(freqs), split = ".", fixed = TRUE), 
    stringsAsFactors = FALSE))[1, ]), bypop = TRUE, overall = TRUE)
}
\arguments{
  \item{freqs}{
A data frame of allele frequencies, such as those output by \code{\link{simpleFreq}} and 
\code{\link{deSilvaFreq}}.
}
  \item{pops}{
An optional characer vector containing names of populations to include.
}
  \item{loci}{
An optional character vector containing names of loci to include.
}
  \item{bypop}{
If \code{TRUE}, PIC will be estimated separately for each population.
}
  \item{overall}{
If \code{TRUE}, mean allele frequencies will be estimated across all populations (weighted
by population size) and used to estimate overall PIC values for each locus.
}
}
\details{
PIC is estimated as:

\deqn{1 - (\sum_{i = 1}^{n} p_{i}^{2}) - \sum_{i = 1}^{n - 1} \sum_{j = i + 1}^{n} 2p_{i}^{2}p_{j}^{2}}{
      1 - (\sum {i = 1} to {n} p_i^2) - \sum {i = 1} to {n - 1} \sum {j = i + 1} to {n} 2 * p_i^2 * p_j^2}

according to Botstein et al. (1980), where \eqn{p_i} and \eqn{p_j} are allele frequencies at alleles 
\emph{i} and \emph{j}, respectively, and \eqn{n} is the number of alleles.

The higher this value is, the more useful a marker is for distinguishing individuals and understanding
relationships among them.
}
\value{
A matrix, with loci in columns, and populations and/or \dQuote{Overall} in rows.  Each element of
the matrix contains a PIC value.
}
\references{
Botstein, M., White, R. L., Skolnick, M. and Davis, R. W. (1980) Construction of 
a genetic linkage map in man using restriction fragment length polymorphisms.  
\emph{American Journal of Human Genetics} \bold{32}, 314--331.
}
\author{
Lindsay V. Clark
}
\seealso{
\code{\link{alleleDiversity}}
}
\examples{
# generate allele frequencies for this example
myfreq <- data.frame(row.names = c("pop1", "pop2"), 
                     Genomes = c(20,30),
                     loc1.124 = c(0.1, 0.25),
                     loc1.126 = c(0.2, 0),
                     loc1.128 = c(0.05, 0.4),
                     loc1.130 = c(0.3, 0.1),
                     loc1.132 = c(0.1, 0.1),
                     loc1.134 = c(0.25, 0.15),
                     loc2.150 = c(0.4, 0.5),
                     loc2.155 = c(0.3, 0.2),
                     loc2.160 = c(0.3, 0.3))

# estimate PIC
PIC(myfreq)

}

\keyword{ arith }
