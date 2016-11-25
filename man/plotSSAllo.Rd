\name{plotSSAllo}
\alias{plotSSAllo}
\alias{plotParamHeatmap}
\alias{processDatasetAllo}
\title{
Perform Allele Assignments across Entire Dataset
}
\description{
\code{processDatasetAllo} runs \code{\link{alleleCorrelations}} on every locus in 
a \code{"\linkS4class{genambig}"} object, then runs \code{\link{testAlGroups}} on every locus
using several user-specified parameter sets.  It chooses a single best set of allele assignments
for each locus, and produces plots to help the user evaluate assignment quality.
\code{plotSSAllo} assists the user in evaluating the quality of allele assignments by plotting
the results of K-means clustering.  \code{plotParamHeatmap} assists the user in choosing the
best parameter set for \code{testAlGroups} for each locus.
}
\usage{
plotSSAllo(AlCorrArray)
plotParamHeatmap(propMat, popname = "AllInd", col = grey.colors(12)[12:1], main = "")
processDatasetAllo(object, samples = Samples(object), loci = Loci(object),
                   n.subgen = 2, SGploidy = 2, n.start = 50, alpha = 0.05,
                   parameters = data.frame(tolerance     = c(0.05, 0.05,  0.05, 0.05),
                                           swap          = c(TRUE, FALSE, TRUE, FALSE),
                                           null.weight   = c(0.5,  0.5,   0,    0)),
                   plotsfile = "alleleAssignmentPlots.pdf", usePops = FALSE, ...)
}

\arguments{
  \item{AlCorrArray}{
A two-dimensional list, where each item in the list is the output of \code{alleleCorrelations}.
The first dimension represents loci, and the second dimension represents populations.  Both 
dimensions are named.  This is the \code{$AlCorrArray} output of \code{processDatasetAllo}.
}
\item{propMat}{A two-dimensional array, with loci in the first dimension and parameter sets in the
second dimension, indicating the proportion of alleles that were found to be homoplasious by \code{testAlGroups}
or the proportion of genotypes that could not be recoded using a given set of allele assignments.
This can be the \code{$propHomoplasious} output of \code{processDatasetAllo}, indexed by a single population.  If a 
three-dimensional array is provided, it will be indexed in the second dimension by \code{popname}.  The 
\code{$propHomoplMerged} or \code{$missRate} output of \code{processDatasetAllo} may also be passed to this argument.}
\item{popname}{The name of the population corresponding to the data in \code{propMat}.}
\item{col}{The color scale for representing the proportion of loci that are homoplasious or the proportion of 
genotypes that are missing.}
\item{main}{A title for the plot.}
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
populations are analyzed separately with \code{alleleCorrelations} and \code{testAlGroups}, before merging the 
results with \code{mergeAlleleAssignments}.}
\item{...}{Additional parameters to pass to \code{testAlGroups} for adjusting the simulated annealing algorithm.}
}
\details{
\code{plotSSAllo} produces a plot of loci by population, with the sums-of-squares ratio on the x-axis and the evenness of allele distribution
on the y-axis (see Value).  Locus names are written directly on the plot.  If there are multiple population names, locus names are colored
by population, and a legend is provided for colors.  Loci with high-quality allele clustering are expected to be in the upper-right
quadrant of the plot.  If locus names are in italics, it indicates that positive correlations were found between some alleles, indicating
population structure or scoring error that could interfere with assignment quality.

\code{plotParamHeatmap} produces an image to indicate the proportion of alleles found to be homoplasious, or the proportion of genotypes that
could not be unambiguously recoded using allele assignments, for each locus and
parameter set for a given population (when looking at homoplasy) or merged across populations (for homoplasy or the proportion of non-recodeable
genotypes).  Darker colors indicate more homoplasy or more genotypes that could not be recoded.  
The word \dQuote{best} indicates, for each
locus, the parameter set that found the least homoplasy or smallest number of non-recodeable genotypes.

By default, \code{processDatasetAllo} generates a PDF file containing output from \code{plotSSAllo} and \code{plotParamHeatmap},
as well as heatmaps of the \code{$heatmap.dist} output of \code{alleleCorrelations} for each locus and population.
Heatmaps are not plotted for loci where an allele is present in all individuals.  \code{processDatasetAllo} also
generates a list of R objects containing allele assignments under different parameters, as well as statistics for evaluating 
clustering quality and choosing the optimal parameter sets, as described below.
}
\value{
\code{plotSSAllo} draws a plot and invisibly returns a list:
\item{ssratio }{A two-dimensional array with loci in the first dimension and populations in the second dimension.
Each value is the sums-of-squares between isoloci divided by the total sums-of-squares, as output by K-means
clustering.  If K-means clustering was not performed, the value is zero.}
\item{evenness }{An array of the same dimensions as \code{$ssratio}, containing values to indicate how evenly
alleles are distributed among isoloci as determined by K-means clustering.  This is:
\deqn{1 - \sum_{1}^{i}{(\frac{a_{i}}{A})^2}}{1 - \sum_{1}^{i}{(a_i/A)^2}}
where \eqn{i} is the number of isoloci, \eqn{a_i} is the number of alleles for a given isolocus, and \eqn{A}
is the total number of alleles for the locus.}
\item{max.evenness }{The maximum possible value for \code{$evenness}, given the number of isoloci.}
\item{min.evenness }{The minimum possible value for \code{$evenness}, given the number of isoloci and alleles.}
\item{posCor }{An array of the same dimensions as \code{$ssratio}, containing \code{TRUE} if there were any
positive correlations between alleles, and \code{FALSE} if not.}

\code{processDatasetAllo } returns a list:
\item{AlCorrArray }{A two-dimensional list with loci in the first dimension and populations in the second
dimension, giving the results of \code{alleleCorrelations}.}
\item{TAGarray }{A three-dimensional list with loci in the first dimension, populations in the second dimension,
and parameter sets in the third dimension, giving the results of \code{testAlGroups}.}
\item{plotSS }{The output of \code{plotSSAllo}.}
\item{propHomoplasious }{A three-dimensional array, with the same dimensions as \code{$TAGarray}, indicating 
the proportion of alleles that were found to be homoplasious for each locus, population, and parameter set.}
\item{mergedAssignments }{A two-dimensional list, with loci in the first dimension and parameter sets in the
second dimension, containing allele assignments merged across populations.  This is the output of 
\code{\link{mergeAlleleAssignments}}.}
\item{propHomoplMerged }{A two-dimensional array, of the same dimensions as \code{$mergedAssignments}, indicating
the proportion of alleles that were homoplasious, for each locus and parameter set, for allele assignments that
were merged across populations.}
\item{missRate }{A matrix with the same dimensions as \code{$mergedAssignments} indicating the proportion of 
non-missing genotypes from the original dataset that cannot be unambiguously recoded, without invoking
aneuploidy, using the merged allele assignments from each parameter set for each locus.}
\item{bestAssign }{A one-dimensional list with a single best set of allele assignments, from \code{$mergedAssignments}, 
for each locus.  The best set of assignments is chosen using \code{$missRate}, then in the case of a tie using 
\code{$propHomoplMerged}, then in the case of a tie using the parameter set that was listed first.}

\code{plotParamHeatmap} draws a plot and does not return anything.
}
\references{
Clark, L. V. and Drauch Schreier, A. (2016) Resolving microsatellite genotype ambiguity in
populations of allopolyploid and diploidized autopolyploid organisms
using negative correlations between allelic variables. \emph{Molecular Ecology Resources}, 
DOI: 10.1111/1755-0998.12639.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{alleleCorrelations}}, \code{\link{recodeAllopoly}}
}
\examples{
# get example dataset
data(AllopolyTutorialData)

# data cleanup
mydata <- deleteSamples(AllopolyTutorialData, c("301", "302", "303"))
PopInfo(mydata) <- rep(1:2, each = 150)
Genotype(mydata, 43, 2) <- Missing(mydata)

# allele assignments
# R is set to 10 here to speed processing for example.  It should typically be left at the default.
myassign <- processDatasetAllo(mydata, loci = c("Loc3", "Loc6"),
                               plotsfile = NULL, usePops = TRUE, R = 10,
                               parameters = data.frame(tolerance = c(0.5, 0.5), 
                               swap = c(TRUE, FALSE),
                               null.weight = c(0.5, 0.5)))

# view best assignments for each locus
myassign$bestAssign
                               
# plot K-means results
plotSSAllo(myassign$AlCorrArray)

# plot proportion of homoplasious alleles
plotParamHeatmap(myassign$propHomoplasious, "Pop1")
plotParamHeatmap(myassign$propHomoplasious, "Pop2")
plotParamHeatmap(myassign$propHomoplMerged, "Merged across populations")

# plot proportion of missing data, after recoding, for each locus and parameter set
plotParamHeatmap(myassign$missRate, main = "Missing data:")
}

\keyword{ hplot }
\keyword{ iteration }
