\name{read.ATetra}
\alias{read.ATetra}
\title{
Read File in ATetra Format
}
\description{
Given a file formatted for the software ATetra, \code{read.ATetra}
produces a \code{genambig} object containing genotypes, population
identities, population names, and a dataset description from the file.
Ploidy in the \code{genambig} object is automatically set to
\code{4}.
}
\usage{
read.ATetra(infile)
}
\arguments{
  \item{infile}{
Character string.  A file path to the file to be read.
}
}
\details{
\code{read.ATetra} reads text files in the exact format specified by
the ATetra documentation.  Note that this format only allows tetraploid
 data and that there can be no missing data.
}
\value{
A \code{genambig} object as described above.
}
\references{
  \url{http://www.vub.ac.be/APNA/ATetra_Manual-1-1.pdf}

  van Puyvelde, K., van Geert, A. and Triest, L.  (2010)  ATETRA, a new
  software program to analyze tetraploid microsatellite data: comparison
  with TETRA and TETRASAT.  \emph{Molecular Ecology Resources}
  \bold{10}, 331--334.
}
\author{
Lindsay V. Clark
}

\seealso{
  \code{\link{write.ATetra}}, \code{\link{read.Tetrasat}},
  \code{\link{read.GeneMapper}},
\code{\link{read.Structure}}, \code{\link{read.GenoDive}},
\code{\link{read.SPAGeDi}}, \code{\link{read.POPDIST}},
\code{\link{read.STRand}}
}
\examples{
# create a file to be read
# (this would normally be done in a text editor or with ATetra's Excel template)
cat("TIT,Sample Rubus Data for ATetra", "LOC,1,CBA15",
"POP,1,1,Commonwealth", "IND,1,1,1,CMW1,197,208,211,213",
"IND,1,1,2,CMW2,197,207,211,212", "IND,1,1,3,CMW3,197,208,212,219",
"IND,1,1,4,CMW4,197,208,212,219", "IND,1,1,5,CMW5,197,208,211,212",
"POP,1,2,Fall Creek Lake", "IND,1,2,6,FCR4,197,207,211,212",
"IND,1,2,7,FCR7,197,208,212,218", "IND,1,2,8,FCR14,197,207,212,218",
"IND,1,2,9,FCR15,197,208,211,212", "IND,1,2,10,FCR16,197,208,211,212",
"IND,1,2,11,FCR17,197,207,212,218","LOC,2,CBA23","POP,2,1,Commonwealth",
"IND,2,1,1,CMW1,98,100,106,125","IND,2,1,2,CMW2,98,125,,",
"IND,2,1,3,CMW3,98,126,,","IND,2,1,4,CMW4,98,106,119,127",
"IND,2,1,5,CMW5,98,106,125,","POP,2,2,Fall Creek Lake",
"IND,2,2,6,FCR4,98,125,,","IND,2,2,7,FCR7,98,106,126,",
"IND,2,2,8,FCR14,98,127,,","IND,2,2,9,FCR15,98,108,117,",
"IND,2,2,10,FCR16,98,125,,","IND,2,2,11,FCR17,98,126,,","END",
file = "atetraexample.txt", sep = "\n")

# Read the file and examine the data
exampledata <- read.ATetra("atetraexample.txt")
summary(exampledata)
PopNames(exampledata)
viewGenotypes(exampledata)
}
\keyword{file}
