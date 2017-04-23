\name{genIndex}
\alias{genIndex}
\alias{genIndex,genambig-method}
\alias{genIndex,array-method}
\title{
Find All Unique Genotypes for a Locus
}
\description{
This function will return all unique genotypes for a given locus (ignoring allele
order, but taking copy number into account) and return those genotypes as well 
as an index indicating which genotype(s) each individual has.  This is a generic function
with methods for \code{"\linkS4class{genambig}"} objects and for arrays.  The array 
method is primarily intended for internal use with \code{\link{meandistance.matrix2}},
processing the output of \code{\link{genotypeProbs}}.
}
\usage{
genIndex(object, locus)
}
\arguments{
  \item{object}{
Typically, a \code{"genambig"} object.  A two-dimentional list (array) can also be
used here, where samples are in the first dimension and loci in the second dimension
and each element of the list is output from \code{genotypeProbs}.
}
  \item{locus}{
A character string or integer indicating which locus to process.
}
}
\value{
A list with two elements:

\item{uniquegen }{A list, where each element in the list is a vector indicating a
unique genotype that was found.}
\item{genindex }{For \code{"genambig"} objects, an integer vector, with one value per sample.
This is the index of that sample's genotype in \code{uniquegen}.  For arrays, a list with
one element per sample.  Each element is a vector of indices of that sample's possible 
genotypes in \code{uniquegen}, in the same order as in the \code{genotypeProbs} output.}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{meandistance.matrix}} uses the \code{"genambig"} method internally.

\code{\link{.unal1loc}}, \code{\link{assignClones}}
}
\examples{
data(simgen)
genIndex(simgen, 1)
}
\keyword{ methods }
\keyword{ manip }
