\name{gendata.to.genind}
\alias{gendata.to.genind}
\title{
Convert Data to genind Format
}
\description{
This is a function for exporting data to the package \pkg{adegenet}.
}
\usage{
gendata.to.genind(object, samples = Samples(object), loci = Loci(object))
}
\arguments{
  \item{object}{
A \code{"genambig"} or (preferably) a \code{"genbinary"} object.  
}
  \item{samples}{
A character vector indicating the samples to include in the output.
}
  \item{loci}{
A character vector indicating the loci to include in the output.
}
}
\details{
\code{gendata.to.genind} converts a \code{"genambig"} or
\code{"genbinary"} object to a \code{"genind"} object using the package
\pkg{adegenet}.  The data being converted must all be of a single
ploidy.  Ploidy and population information are carried over to the new
object.  Data will be coded as presence/absence in the new object.  The
locus names in the new object are locus and allele names
seperated by a hyphen.

\pkg{adegenet} must be installed in order to use this function.
}
\value{
A genetic dataset in the \code{"genind"} class, ready for use in
\pkg{adegenet}. 
}
\references{
  \url{http://adegenet.r-forge.r-project.org/}
  
Jombart, T.  (2008)  adegenet: a R package for the multivariate analysis
of genetic markers.  \emph{Bioinformatics} \bold{24}, 1403--1405.
}
\author{
Lindsay V. Clark
}


\seealso{
\code{\link{freq.to.genpop}}
}
\examples{
require(adegenet)

# create a "genambig" object
mydata <- new("genambig", samples=c("a","b","c","d"), loci=c("e","f"))
PopNames(mydata) <- c("G","H")
PopInfo(mydata) <- c(1,1,2,2)
mydata <- reformatPloidies(mydata, output="one")
Ploidies(mydata) <- 3
Genotypes(mydata, loci="e") <- list(c(100),c(100,102),
                                    c(98,102,104),c(102,106))
Genotypes(mydata, loci="f") <- list(c(200,202,204),Missing(mydata),
                                    c(210,212),c(204,210,214))

# convert to "genind"
mydata2 <- gendata.to.genind(mydata)
mydata2@tab
locNames(mydata2)
indNames(mydata2)
popNames(mydata2)
pop(mydata2)
}
\keyword{ manip }

