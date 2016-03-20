\name{plotSSAllo}
\alias{plotSSAllo}
\alias{processDatasetAllo}
\title{
Perform Allele Asssignments across Entire Dataset
}
\description{
\code{processDatasetAllo} runs \code{\link{alleleCorrelations}} on every locus in 
a \code{"\linkS4class{genambig}"} object, then runs \code{\link{testAlGroups}} on every locus
using several user-specified parameter sets.  % additional functions; picking parameters, plots
\code{plotSSAllo} assists the user in evaluating the quality of allele assignments by plotting
the results of K-means clustering.
}
\usage{
plotSSAllo(AlCorrArray)
processDatasetAllo(object, samples = Samples(object), loci = Loci(object),
                   n.subgen = 2, SGploidy = 2, n.start = 50, alpha = 0.05,
                   parameters = data.frame(tolerance     = c(0.05, 0.05, 0.05),
                                           rare.al.check = c(0.2,  0,    0.2),
                                           null.weight   = c(0.5,  0.5,  0)),
                   plotsfile = "alleleAssignmentPlots.pdf", usePops = FALSE)
}

\arguments{
  \item{AlCorrArray}{
A two-dimensional list, where each item in the list is the output of \code{alleleCorrelations}.
The first dimension represents loci, and the second dimension represents populations.  Both 
dimensions are named.  This is the \code{$AlCorrArray} output of \code{processDatasetAllo}.
}
\item{object}{A \code{"genambig"} object.}
\item{samples}{An optional character vector indicating which samples to include in analysis.}
\item{loci}{An optional character vector indicating which loci to include in analysis.}
\item{n.subgen}{The number of isoloci into which each locus should be split.  Passed directly to 
\code{alleleCorrelations}.}
\item{SGploidy}{The ploidy of each isolocus.  Passed directly to \code{testAlGroups}.}
\item{n.start}{Passed directly to the \code{nstart} argument of \code{kmeans}.  See \code{\link{alleleCorrelations}}.}
\item{alpha}{The significance threshold for determining whether two alleles are significantly correlated.  Used
primarily for identifying potentially problematic positive correlations.  Passed directly to \code{alleleCorrelations}.}
\item{parameters}{Data frame indicating parameter sets to pass to \code{testAlGroups}.  Each row is one set of parameters.}
\item{plotsfile}{A PDF output file name for drawing plots to help assess assignment quality.  Can be \code{NULL} if no
plots are desired.}
\item{usePops}{If \code{TRUE}, population assignments are taken from the \code{PopInfo} slot of \code{object}, and 
populations are analyzed separately.}
}
\details{
\code{plotSSAllo} produces a plot of loci by population, with the sums-of-squares ratio on the x-axis and the evenness of allele distribution
on the y-axis (see Value).  Locus names are written directly on the plot.  If multiple population names, locus names are colored
by population, and a legend is provided for colors.  Loci with high-quality allele clustering are expected to be in the upper-right
quadrant of the plot.
}
\value{
\code{plotSSAllo} draws a plot and invisibly returns a list:
\item{ssratio }{A two-dimensional array with loci in the first dimension and populations in the second dimension.
Each value is the sums-of-squares between isoloci divided by the total sums-of-squares, as output by K-means
clustering.}
\item{evenness }{An array of the same dimensions as \code{$ssratio}, containing values to indicate how evenly
alleles are distributed among isoloci as determined by K-means clustering.  This is:
\deqn{1 - \sum_{1}^{i}{(frac{a_{i}}{A})^2}}{1 - \sum_{1}^{i}{(a_i/A)^2}}
where \eqn{i} is the number of isoloci, \eqn{a_i} is the number of alleles for a given isolocus, and \eqn{A}
is the total number of alleles for the locus.}
\item{max.evenness }{The maximum possible value for \code{$evenness}, given the number of isoloci.}
\item{min.evenness }{The minimum possible value for \code{$evenness}, given the number of isoloci and alleles.}

\code{processDatasetAllo } returns a list:
\item{AlCorrArray }{A two-dimensional list with loci in the first dimension and populations in the second
dimension, giving the results of \code{alleleCorrelations}.}
\item{TAGarray }{A three-dimensiona list with loci in the first dimension, populations in the second dimension,
and parameter sets in the third dimension, giving the results of \code{testAlGroups}.}
% additional items as the function is added to.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# get example dataset
dataset(AllopolyTutorialData)

# data cleanup
mydata <- deleteSamples(AllopolyTutorialData, c("301", "302", "303"))
PopInfo(mydata) <- rep(1:2, each = 150)
Genotype(mydata, 43, 2) <- Missing(mydata)

# allele assignments
myassign <- processDatasetAllo(mydata, loci = c("Loc3", "Loc6"),
                               plotsfile = NULL, usePops = TRUE)
# plot K-means results
plotSSAllo(myassign$AlCorrArray)

}

\keyword{ hplot }
\keyword{ iteration }
