\name{read.Tetrasat}
\alias{read.Tetrasat}
\title{
Read Data from a TETRASAT Input File
}
\description{
Given a file containing genotypes in the TETRASAT format,
\code{read.Tetrasat} produces a \code{genambig} object containing
genotypes and population identities from the file.
}
\usage{
read.Tetrasat(infile)
}
\arguments{
  \item{infile}{
A character string of the file path to be read.
}
}
\details{
\code{read.Tetrasat} reads text files that are in the exact format
specified by the software TETRASAT  and TETRA (see references for more
information).  This is similar to the file format for GenePop but allows
for up to four alleles per locus.  All alleles must be coded by two
digits.  Another difference between the TETRASAT and GenePop formats is
that in TETRASAT the sample name and genotypes are not separated by a
comma, because the columns of data have fixed widths.

Since TETRASAT files also contain information about which samples belong
to which populations, this information is put into the \code{PopInfo}
slot of the \code{genambig} object.  Population names are not taken from
the file.  The \code{Ploidies} slot is
filled with the number 4 (using the \code{"\linkS4class{ploidyone}"}
class), because all individuals should be tetraploid.
The first line of the file is put into the \code{Description} slot.

Locus names should not contain the letters "pop", uppercase or lowercase,
adjacent to each other.
}
\value{
A \code{genambig} object containing data from the file.
}
\references{
  Markwith, S. H., Stewart, D. J. and Dyer, J. L.  (2006)
  TETRASAT: a program for the population analysis of
  allotetraploid microsatellite data. \emph{Molecular Ecology Notes}
  \bold{6}, 586-589.

%  \url{http://ecology.bnu.edu.cn/zhangdy/TETRA/TETRA.htm} % zhangdy/index.php still up, but TETRA link broken July 2016

  Liao, W. J., Zhu, B. R., Zeng, Y. F. and Zhang, D. Y.  (2008)
  TETRA: an improved program for population genetic analysis of
  allotetraploid microsatellite data.  \emph{Molecular Ecology Resources}
  \bold{8}, 1260--1262.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{read.GeneMapper}}, \code{\link{write.Tetrasat}},
\code{\link{read.ATetra}}, \code{\link{read.GenoDive}},
\code{\link{read.Structure}},
\code{\link{read.SPAGeDi}}, \code{\link{read.POPDIST}},
\code{\link{read.STRand}}
}
\examples{
# example with defined data:
myfile <- tempfile()
cat("Sample Data", "A1_Gtype", "A10_Gtype", "B1_Gtype", "D7_Gtype",
"D9_Gtype", "D12_Gtype", "Pop",
"BCRHE 1             0406     04040404 0208     02020202 03030303 0710",
"BCRHE 10            0406     04040404 07070707 02020202 0304     0710",
"BCRHE 2             04040404 04040404 0708     02020202 010305   0710",
"BCRHE 3             04040404 04040404 02020202 0203     03030303 0809",
"BCRHE 4             04040404 04040404 0608     0203     03030303 070910",
"BCRHE 5             04040404 04040404 0208     02020202 03030303 050710",
"BCRHE 6             0304     04040404 0207     02020202 03030303 07070707",
"BCRHE 7             0406     04040404 0708     02020202 03030303 07070707",
"BCRHE 8             0304     04040404 0203     0203     03030303 0709",
"BCRHE 9             0406     04040404 0708     02020202 03030303 0710",
"Pop",
"BR 1                0406     04040404 05050505 02020202 03030303 1012",
"BR 10               030406   04040404 0607     02020202 03030303 1011",
"BR 2                030406   04040404 07070707 02020202 03030303 09090909",
"BR 3                010304   04040404 07070707 02020202 03030303 09090909",
"BR 4                030406   04040404 07070707 0203     03030303 10101010",
"BR 5                030406   04040404 07070707 02020202 03030303 10101010",
"BR 6                0406     04040404 0507     0203     03030303 10101010",
"BR 7                0304     04040404 0809     02020202 03030303 070910",
"BR 8                030406   04040404 07070707 02020202 03030303 070910",
"BR 9                0406     04040404 07070707 02020202 03030303 07070707",
sep="\n", file=myfile)
mydata2 <- read.Tetrasat(myfile)

summary(mydata2)
viewGenotypes(mydata2, loci="B1_Gtype")
}
\keyword{file}
