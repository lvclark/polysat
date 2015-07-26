\name{read.STRand}
\alias{read.STRand}
\title{
Read Genotypes Produced by STRand Software
}
\description{
This function reads in data in a format derived from the \dQuote{BTH}
format for exporting genotypes from the allele calling software STRand.
}
\usage{
read.STRand(file, sep = "\t", popInSam = TRUE)
}
\arguments{
  \item{file}{
A text string indicating the file to read.
}
  \item{sep}{
Field delimiter for the file.  Tab by default.
}
  \item{popInSam}{
Boolean.  If \code{TRUE}, fields from the \dQuote{Pop} and \dQuote{Ind}
columns will be concatenated to create a sample name.  If \code{FALSE},
only the \dQuote{Ind} column will be used for sample names.
}
}
\details{
This function does not read the files directly produced from STRand, but
requires some simple clean-up in spreadsheet software.  The BTH format
in STRand produces two columns per locus.  One of these columns should
be deleted so that there is just one column per locus.  Loci names
should remain in the column headers.  The column containing sample names
should be deleted or renamed \dQuote{Ind}.  A
\dQuote{Pop} column will need to be added, containing population names.
An \dQuote{Ind} column is also necessary, containing either full sample
names or a sample suffix to be concatenated with the population name
(see \code{popInSam} argument).

STRand adds an asterisk to the end of any genotype with more than two
alleles.  \code{read.STRand} will automatically strip this asterisk out
of the genotype.

Missing data is indicated by a zero in the file.
}
\value{
A \code{"genambig"} object containing genotypes, locus and sample names,
population names, and population identities from the file.
}
\references{
  \url{http://www.vgl.ucdavis.edu/informatics/strand.php/}

  Toonen, R. J. and Hughes, S. (2001) Increased Throughput for Fragment
  Analysis on ABI Prism 377 Automated Sequencer Using a Membrane Comb
  and STRand Software.  \emph{Biotechniques} \bold{31}, 1320--1324.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.table}}, \code{\link{read.GeneMapper}},
\code{\link{read.GenoDive}}, \code{\link{read.Structure}},
\code{\link{read.ATetra}}, \code{\link{read.Tetrasat}},
\code{\link{read.POPDIST}}, \code{\link{read.SPAGeDi}}
}
\examples{
# generate file to read
strtemp <- data.frame(Pop=c("P1","P1","P2","P2"),
                      Ind=c("a","b","a","b"),
                      LocD=c("0","172/174","170/172/178*","172/176"),
                      LocG=c("130/136/138/142*","132/136","138","132/140/144*"))
write.table(strtemp, file="strtemp.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

# read the file
mydata <- read.STRand("strtemp.txt")
viewGenotypes(mydata)
PopNames(mydata)

# alternative example with popInSam=FALSE
strtemp$Ind <- c("OH1","OH5","MT4","MT7")
write.table(strtemp, file="strtemp.txt", sep="\t",
            row.names=FALSE, quote=FALSE)
mydata <- read.STRand("strtemp.txt", popInSam=FALSE)
Samples(mydata)
PopNames(mydata)
}
\keyword{ file }
