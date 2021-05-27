\name{write.freq.SPAGeDi}
\alias{write.freq.SPAGeDi}
\title{
Create a File of Allele Frequencies for SPAGeDi
}
\description{
A table of allele frequencies such as that produced by \code{simpleFreq}
or \code{deSilvaFreq} is used to calculate average allele frequencies
for the entire dataset.  These are then written in a format that can be
read by the software SPAGeDi.
}
\usage{
write.freq.SPAGeDi(freqs, usatnts, file = "", digits = 2,
                   pops = row.names(freqs),
  loci = unique(as.matrix(as.data.frame(strsplit(names(freqs), split =
  ".", fixed = TRUE), stringsAsFactors = FALSE))[1, ]))
}
\arguments{
  \item{freqs}{
A data frame of population sizes and allele frequencies, such as that
produced by \code{simpleFreq} or \code{deSilvaFreq}.  Populations are in
rows, and alleles are in columns.  A column is needed containing
population sizes in number of genomes; this may either be a single
column called \dQuote{Genomes} or multiple columns named by the locus
and \dQuote{Genomes}, sepearated by a period.  All other columns contain
allele frequencies.  The column
names for these should be the locus and allele separated by a period.
}
  \item{usatnts}{
An integer vector containing the lengths of the microsatellite repeats
for the loci in the table.  In most cases, if \code{object} is the
\code{"gendata"} object used to generate \code{freqs}, then you should set
\code{usatnts = Usatnts(object)}.  This is needed to convert allele
names in the same way that \code{write.SPAGeDi} converts allele names.
}
  \item{file}{
The name of the file to write.
}
  \item{digits}{
The number of digits to use to represent each allele.  This should be
the same as that used in \code{write.SPAGeDi}, so that allele names are
consistent between the two files.
}
  \item{pops}{
An optional character vector indicating a subset of populations from the
table to use in calculating mean allele frequencies.
}
  \item{loci}{
An optional character vector indicating a subset of loci to write to the
file.
}
}
\details{
For some calculations of inter-individual relatedness and kinship
coefficients, SPAGeDi can read a file of allele frequencies to use in
the calculation.  \code{write.freq.SPAGeDi} puts allele frequencies from
\pkg{polysat} into this format.

A weighted average of allele frequencies is calculated across all
populations (or those specified by \code{pops}).  The average is
weighted by population size as specified in the \dQuote{Genomes} column
of \code{freqs}.

Allele names are converted to match those produced by
\code{write.SPAGeDi}.  Alleles are divided by the numbers in
\code{usatnts} in order to convert fragment length in nucleotides to
repeat numbers.  If necessary, \code{10^(digits-1)} is repeatedly
subtracted from all alleles until they can be represented using the
right number of digits.

The file produced is tab-delimited and contains two columns per locus.
The first column contains the locus name followed by all allele names,
and the second column contains the number of alleles followed by the
allele frequencies.
}
\value{
A file is written but no value is returned.
}
\references{
    \url{https://ebe.ulb.ac.be/ebe/SPAGeDi.html}

  Hardy, O. J. and Vekemans, X. (2002) SPAGeDi: a versatile computer
  program to analyse spatial genetic structure at the individual or
  population levels.  \emph{Molecular Ecology Notes} \bold{2}, 618-620.
}
\author{
Lindsay V. Clark
}
\note{
SPAGeDi can already estimate allele frequencies in a way that is
  identical to that of \code{simpleFreq}.  Therefore, if you have allele
  frequencies produced by \code{simpleFreq}, there is not much sense in
  exporting them to SPAGeDi.  \code{deSilvaFreq}, however, is a more
  advanced and accurate allele frequency estimation than what is
  available in SPAGeDi v1.3.  \code{write.freq.SPAGeDi} exists primarily to
  export allele frequencies from \code{deSilvaFreq}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{write.SPAGeDi}}, \code{\link{deSilvaFreq}}
}
\examples{
\dontrun{
# set up a genambig object to use in this example
mygen <- new("genambig", samples=c(paste("G", 1:30, sep=""),
                                   paste("R", 1:50, sep="")),
             loci=c("afrY", "ggP"))
PopNames(mygen) <- c("G", "R")
PopInfo(mygen) <- c(rep(1, 30), rep(2, 50))
mygen <- reformatPloidies(mygen, output="one")
Ploidies(mygen) <- 4
Usatnts(mygen) <- c(2, 2)

# randomly create genotypes according to pre-set allele frequencies
for(s in Samples(mygen, populations=1)){
    Genotype(mygen, s, "afrY") <-
        unique(sample(c(140, 142, 146, 150, 152), 4, TRUE,
                      c(.30, .12, .26, .08, .24)))
    Genotype(mygen, s, "ggP") <-
        unique(sample(c(210, 214, 218, 220, 222), 4, TRUE,
                      c(.21, .13, .27, .07, .32)))
}
for(s in Samples(mygen, populations=2)){
    Genotype(mygen, s, "afrY") <-
        unique(sample(c(140, 142, 144, 150, 152), 4, TRUE,
                      c(.05, .26, .17, .33, .19)))
    Genotype(mygen, s, "ggP") <-
        unique(sample(c(212, 214, 220, 222, 224), 4, TRUE,
                      c(.14, .04, .36, .20, .26)))
}

# write a SPAGeDi file
write.SPAGeDi(mygen, file="SPAGdataFreqExample.txt")

# calculate allele frequenies
myfreq <- deSilvaFreq(mygen, self = 0.05)

# write allele frequencies file
write.freq.SPAGeDi(myfreq, usatnts=Usatnts(mygen),
file="SPAGfreqExample.txt")
}
\dontshow{
# a shorter example for R CMD check
myfreq2 <- data.frame(row.names=c("G","R"), Genomes=c(120,200),
             abc.130=c(0.13, 0.2), abc.132=c(0.34, 0.5),
             abc.134=c(0.53, 0.3))
myusatnts <- as.integer(2)
names(myusatnts) <- "abc"
write.freq.SPAGeDi(myfreq2, usatnts=myusatnts, file="SPAGfreqExample2.txt")
}

}
\keyword{ file }
\keyword{ array }
