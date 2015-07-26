\name{freq.to.genpop}
\alias{freq.to.genpop}
\title{
Convert Allele Frequencies for Adegenet
}
\description{
Given a data frame of allele frequencies such as that produced by
\code{simpleFreq} or \code{deSilvaFreq}, \code{freq.to.genpop} creates a
data frame of allele counts that can be read by the \code{as.genpop}
function in the package \pkg{adegenet}.
}
\usage{
freq.to.genpop(freqs, pops = row.names(freqs),
               loci =
                unique(as.matrix(as.data.frame(strsplit(names(freqs),
                split = ".", fixed = TRUE), stringsAsFactors = FALSE))[1, ]))
}
\arguments{
  \item{freqs}{
A data frame of allele frequencies.  Row names are population names.
The first column is called \code{"Genomes"} and indicates the size of
each population in terms of number of haploid genomes.  All other column
names are the locus and allele separated by a period.  These columns
contain the frequencies of each allele in each population.  For each
locus and population, all frequencies should total to 1.
}
  \item{pops}{
An optional character vector indicating the names of populations to use.
}
  \item{loci}{
An optional character vector indicating the names of loci to use.
}
}
\details{
  \pkg{adegenet} expects one ploidy for the entire dataset.  Therefore, data
  frames of allele frequencies with multiple \dQuote{Genomes} columns,
  such as those produced when ploidy varies by locus, are not allowed as
  the \code{freqs} argument.
}
\value{
A data frame with row and column names identical to those in
\code{freqs}, minus the \code{"Genomes"} column and any columns for loci
not included in \code{loci}.  Allele frequencies are converted to counts
by multiplying by the values in the \code{"Genomes"} column and rounding
to the nearest integer.
}
\references{
Jombart, T.  (2008)  adegenet: a R package for the multivariate analysis
of genetic markers.  \emph{Bioinformatics} \bold{24}, 1403-1405.
}
\author{
Lindsay V. Clark
}


\seealso{
\code{\link{simpleFreq}}, \code{\link{deSilvaFreq}},
\code{\link{write.freq.SPAGeDi}}, \code{\link{gendata.to.genind}}
}
\examples{
# create a simple allele frequency table
# (usually done with simpleFreq or deSilvaFreq)
myfreq <- data.frame(row.names=c("popA","popB"), Genomes=c(120,100),
                     locG.152=c(0.1,0.4), locG.156=c(0.5, 0.3),
                     locG.160=c(0.4, 0.3), locK.179=c(0.15, 0.25),
                     locK.181=c(0.35, 0.6), locK.183=c(0.5, 0.15))
myfreq

# convert to adegenet format
gpfreq <- freq.to.genpop(myfreq)
gpfreq

\dontrun{
# If you have adegenet installed, you can now make this into a
# genpop object.
require(adegenet)
mygenpop <- genpop(gpfreq, ploidy=as.integer(4), type="codom")

# examine the object that has been created
mygenpop
popNames(mygenpop)
mygenpop@tab
mygenpop@all.names

# Perform a distance calculation with the object
dist.genpop(mygenpop)
}
}
\keyword{ manip }
