\name{AllopolyTutorialData}
\alias{AllopolyTutorialData}
\docType{data}
\title{
Simulated Allotetraploid Data
}
\description{
This is a simulated microsatellite dataset for seven loci and 303 individuals.  
It is intended to be used
with the tutorial \dQuote{Assigning alleles to isoloci in \pkg{polysat}}.
}
\usage{data("AllopolyTutorialData")}
\format{
  A \code{"\linkS4class{genambig}"} object.
}

\examples{
data(AllopolyTutorialData)
summary(AllopolyTutorialData)
viewGenotypes(AllopolyTutorialData, samples=1:10, loci=1)
}
\keyword{datasets}
