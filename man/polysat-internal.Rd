\name{Internal Functions}
\alias{G}
\alias{INDEXG}
\alias{GENLIST}
\alias{RANMUL}
\alias{SELFMAT}
\alias{.unal1loc}
\alias{fixloci}
\title{Internal Functions in polysat}
\description{
  The internal functions \code{G}, \code{INDEXG}, \code{GENLIST},
  \code{RANMUL}, and \code{SELFMAT}
  are used for calculating
  genotype probabilities under partial selfing.  The internal function
  \code{.unal1loc} finds all unique alleles at a single locus.  The internal
  function \code{fixloci} converts locus names to a format that will be
  compatible as column headers for allele frequency tables.
}
\usage{
G(q, n)
INDEXG(ag1, na1, m2)
GENLIST(ng, na1, m2)
RANMUL(ng, na1, ag, m2)
SELFMAT(ng, na1, ag, m2)
.unal1loc(object, samples, locus)
fixloci(loci, warn = TRUE)
}
\arguments{
  \item{q}{Integer.}
  \item{n}{Integer.}
  \item{ag1}{A vector representing an unambiguous genotype.}
  \item{na1}{Integer. The number of alleles, including a null.}
  \item{m2}{Integer.  The ploidy.}
  \item{ng}{Integer.  The number of genotypes.}
  \item{ag}{An array of genotypes such as that produced by
    \code{.genlist}.}
  \item{object}{A \code{"genambig"} object.}
  \item{samples}{Optional, a numeric or character vector indicating
    which samples to use.}
  \item{locus}{A character string or number indicating which locus to
    use.}
  \item{loci}{A character vector of locus names.}
  \item{warn}{Boolean indicating whether a warning should be issued if 
  locus names are changed.}
}
\value{
  \code{G} returns

  \deqn{\frac{(n+q)!}{(q+1)! * (n-1)!}}{(n+q)!/((q+1)! * (n-1)!)}

  \code{INDEXG} returns an integer indicating the row containing a
  particular genotype in the matrix produced by \code{GENLIST}.

  \code{GENLIST} returns an array with dimensions \code{ng, m2},
  containing all possible unambiguous genotypes, one in each row.  The null
  allele is the highest-numbered allele.

  \code{RANMUL} returns a list.  The first item is a vector of
  polynomial coefficients for calculating genotype frequencies under
  random mating.  The second is an array
  showing how many copies of each allele each genotype has.

  \code{SELFMAT} returns the selfing matrix.  Parental genotypes are
  represented in rows, and offspring genotypes in columns.  The numbers
  indicate relative amounts of offspring genotypes produced when the
  parental genotypes are self-fertilized.

  \code{.unal1loc} returns a vector containing all unique alleles, not
  including \code{Missing(object)}.

  \code{fixloci} returns a character vector of corrected locus names.
}

\references{
  De Silva, H. N., Hall, A. J., Rikkerink, E., and Fraser, L. G. (2005)
  Estimation of allele frequencies in polyploids under
certain patterns of inheritance.  \emph{Heredity} \bold{95}, 327--334
}
\seealso{
  \code{\link{deSilvaFreq}}, \code{\link{meandistance.matrix2}},
  \code{\link{genotypeProbs}}, \code{\link{genambig.to.genbinary}},
  \code{\link{alleleDiversity}}
}
\examples{
# Calculation of genotype probabilities in a tetraploid with four
# alleles plus a null, and a selfing rate of 0.5.  This is a translation
# of code in the supplementary material of De Silva et al. (2005).
m2 <- 4
m <- m2/2
na1 <- 5
self <- 0.5
ng <- na1
for(j in 2:m2){
    ng <- ng*(na1+j-1)/j
}
ag <- polysat:::GENLIST(ng, na1, m2)
temp <- polysat:::RANMUL(ng, na1, ag, m2)
rmul <- temp[[1]]
arep <- temp[[2]]
rm(temp)
smat <- polysat:::SELFMAT(ng, na1, ag, m2)
smatdiv <- (polysat:::G(m-1,m+1))^2
p1 <- c(0.1, 0.4, 0.2, 0.2, 0.1) # allele frequencies

# GPROBS subroutine
rvec <- rep(0,ng)
for(g in 1:ng){
    rvec[g] <- rmul[g]
    for(j in 1:m2){
        rvec[g] <- rvec[g]*p1[ag[g,j]]
    }
}
id <- diag(nrow=ng)
smatt <- smat/smatdiv
s3 <- id - self * smatt
s3inv <- solve(s3)
gprob <- (1-self) * s3inv \%*\% rvec
# gprob is a vector of probabilities of the seventy genotypes.
}
\author{Lindsay V. Clark}
\keyword{distribution}
