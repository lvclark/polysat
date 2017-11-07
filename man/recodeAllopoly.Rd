\name{recodeAllopoly}
\alias{recodeAllopoly}
\title{
Create a New \code{\linkS4class{genambig}} Dataset with Loci Split into Isoloci
}
\description{
Given a \code{"\linkS4class{genambig}"} object and a list of allele
assignments such as those produced by \code{\link{testAlGroups}} or
\code{\link{catalanAlleles}}, \code{recodeAllopoly} will generate a
new \code{"genambig"} object, with genotypes split according to which
alleles belong to which isoloci.
}
\usage{
recodeAllopoly(object, x, allowAneuploidy = TRUE,
               samples = Samples(object), loci = Loci(object))
}
\arguments{
  \item{object}{
A \code{"genambig"} object containing the dataset that needs to be re-coded.
}
  \item{x}{
A list.  Each item in the list should itself be a list, in the format
output by \code{testAlGroups}, \code{catalanAlleles}, or
\code{\link{mergeAlleleAssignments}}.  Each sub-list
has three items: \code{$locus} is the name of the locus,
\code{$SGploidy} is an integer indicating the ploidy of each subgenome
(\emph{e.g.} \code{2} for an allotetraploid), and \code{$assignments} is
a matrix of ones and zeros indicating which alleles belong to which isoloci.
}
  \item{allowAneuploidy}{
Boolean.  This controls what happens when the function encounters
genotypes that have more alleles than are possible for a given
isolocus.  (For example, the genotype has four alleles, but three belong
to isolocus 1 and one belongs to isolocus 2.)  If \code{TRUE}, the
individual is assumed to be aneuploid at that locus, and its ploidy is
adjusted only for that locus.  If \code{FALSE}, missing data are recorded.
}
  \item{samples}{
An optional character vector indicating which samples to analyze and output.
}
  \item{loci}{
An optional character vector indicating which loci to analyze and output.
}
}
\details{
The same locus may appear more than once in \code{x}, for example if
distinct populations were analyzed separately to produce the
allele assignments.  If this is the case, \code{recodeAllopoly} will
internally use \code{mergeAlleleAssignments} to consolidate items in
\code{x} with the same locus name.  Loci that are in \code{x} but not
\code{object} are ignored with a warning.  Loci that are in
\code{object} but not \code{x} are retained in the output of the
function, but not re-coded.

This function allows homoplasy, and uses process-of-elimination to try
to determine which isoloci the homoplasious alleles belong to.  In cases
where genotypes cannot be determined for certain due to homoplasy,
missing data are inserted.

If a genotype has more alleles than should be possible (\emph{e.g.} five
alleles in an allotetraploid), the genotype is skipped and will be
output as missing data for all corresponding isoloci.
}
\value{
A \code{"genambig"} object, with loci that are in \code{x} split into
the appropriate number of isoloci.
}
\references{
Clark, L. V. and Drauch Schreier, A. (2017) Resolving microsatellite genotype ambiguity in
populations of allopolyploid and diploidized autopolyploid organisms
using negative correlations between alleles. \emph{Molecular Ecology Resources}, 
\bold{17}, 1090--1103.  DOI: 10.1111/1755-0998.12639.
}
\author{
Lindsay V. Clark
}

\examples{
# generate a dataset for this example
testdata <- new("genambig", samples = paste("S", 1:9, sep = ""),
                loci = c("L1", "L2","L3"))
Genotypes(testdata, loci="L1") <-
    list(c(120,124),c(124,126,130),c(120,126),c(126,132,134),
         c(120,124,130,132),c(120,126,130),c(120,132,134),
         c(120,124,126,130),c(120,132,138))
Genotypes(testdata, loci="L2") <-
    list(c(210,219,222,225),c(216,228),c(210,213,219,222),c(213,222,225,228),
         c(210,213,216,219),c(222,228),c(213),c(210,216),c(219,222,228))
Genotypes(testdata, loci="L3") <-
    list(c(155,145,153),c(157,155),c(151,157,159,165),c(147,151),c(149,153,157),
         c(149,157),c(153,159,161),c(163,165),c(147,163,167))
viewGenotypes(testdata)

# generate allele assignments for this example
myAssign <- list(list(locus="L1", SGploidy=2,
             assignments=matrix(c(1,0,0,1,1,1,0,1,1,0,1,1), nrow=2,
                                ncol=6, dimnames=list(NULL,
                                  c("120","124","126","130","132","134")))),
    list(locus="L2", SGploidy=2,
         assignments=matrix(c(1,1,1,1,1,1,1,0,1,0,1,0,0,1), nrow=2, ncol=7,
    dimnames=list(NULL,c("210","213","216","219","222","225","228")))),
    list(locus="L3", SGploidy=2, assignments="No assignment"))
myAssign

# recode the dataset
splitdata <- recodeAllopoly(testdata, myAssign)

# view results
viewGenotypes(splitdata)
Ploidies(splitdata)

}
\keyword{ manip }

