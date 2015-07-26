\name{Accessors}
\alias{Samples}
\alias{Samples<-}
\alias{Loci}
\alias{Loci<-}
\alias{PopInfo}
\alias{PopInfo<-}
\alias{PopNames}
\alias{PopNames<-}
\alias{Ploidies}
\alias{Ploidies<-}
\alias{Usatnts}
\alias{Usatnts<-}
\alias{Description}
\alias{Description<-}
\alias{Missing}
\alias{Missing<-}
\alias{Genotype}
\alias{Genotype<-}
\alias{Genotypes}
\alias{Genotypes<-}
\alias{PopNum}
\alias{PopNum<-}
\alias{Present}
\alias{Present<-}
\alias{Absent}
\alias{Absent<-}
\title{
Accessor and Replacement Functions for "gendata" Objects
}
\description{
The accessor functions return information that is contained, either
directly or indirectly, in the slots of a \code{gendata} object.  The
replacement functions alter information in one or more slots as appropriate.
}
\usage{
Samples(object, populations, ploidies)
Samples(object) <- value
Loci(object, usatnts, ploidies)
Loci(object) <- value
PopInfo(object)
PopInfo(object) <- value
PopNames(object)
PopNames(object) <- value
PopNum(object, popname)
PopNum(object, popname) <- value
Ploidies(object, samples, loci)
Ploidies(object) <- value
Usatnts(object)
Usatnts(object) <- value
Description(object)
Description(object) <- value
Missing(object)
Missing(object) <- value
Present(object)
Present(object) <- value
Absent(object)
Absent(object) <- value
Genotype(object, sample, locus)
Genotype(object, sample, locus) <- value
Genotypes(object, samples = Samples(object), loci = Loci(object))
Genotypes(object, samples = Samples(object), loci = Loci(object)) <- value
}
\arguments{
  \item{object}{
An object of the class \code{gendata} or one of its subclasses.
}
  \item{populations}{
A character or numeric vector indicating from which populations to return
samples. (optional)
}
  \item{ploidies}{
A numeric vector indicating ploidies, if only samples or loci with a certain
ploidy should be returned.  (optional)
}
\item{sample}{A character string or number indicating the name or number
  of the sample whose genotype should be returned.}
\item{locus}{A character string or number indicating the name or number
  of the locus whose genotype should be returned.}
\item{samples}{A character or numeric vector indicating samples for
  which to return genotypes or ploidies. (optional)}
\item{loci}{A character or numeric vector indicating loci for which to
  return genotypes or ploidies.  (optional)}
\item{usatnts}{A numeric vector indicating microsatellite repeat
  lengths, where only loci of those repeat lengths should be returned.
  (optional)}
\item{popname}{Chacter string or vector.  The name(s) of the
  population(s) for which to retrieve or replace the corresponding
  \code{PopInfo} number(s).  The replacement function should only be
  used for one population at a time.}
\item{value}{\itemize{
  \item For \code{Samples}: a character vector of sample names.
  \item For \code{Loci}: a character vector of locus names.
  \item For \code{PopInfo}: A numeric vector (integer or can be coerced to
  integer) indicating the population identities of samples.
  \item For \code{PopNames}: A character vector indicating the names of
  populations.
  \item For \code{PopNum}: A number (integer or can be coerced to integer)
  that should be the new population number associated with
  \code{popname}.
  \item For \code{Ploidies}: A numeric vector or matrix (integer or can be coerced to
  integer) indicating the ploidy of each sample, locus, and/or the
  dataset.  See \code{\link{reformatPloidies}} and
  \code{"\linkS4class{ploidysuper}"}.
  \item For \code{Usatnts}: A numeric vector (integer or can be coerced to
  integer) indicating the repeat type of each microsatellite locus.
  Dinucleotide repeats should be represented with \code{2},
  trinucleotide repeat with \code{3}, and so on.  If the alleles for a
  given locus are already stored in terms of repeat number rather than
  fragment length in nucleotides, the \code{Usatnts} value for that
  locus should be \code{1}.
  \item For \code{Description}: A character string or character vector
  describing the dataset.
  \item For \code{Missing}: A symbol (usually an integer) to be used to
  indicate missing data.
  \item For \code{Present}: A symbol (usually an integer) to be used to
  indicate the presence of an allele.
  \item For \code{Absent}: A symbol (usually an integer) to be used to
  indicate the absence of an allele.
  \item For \code{Genotype}: a vector of alleles, if the object is of class
  \code{genambig}.
  \item For \code{Genotypes}: A list of vectors (genotypes), of the same
  dimensionality as \code{c(length(samples), length(loci))}, if the
  object is of class \code{genambig}.  If the object is of class
  \code{genbinary}, \code{value} should be a matrix, with column names
  of the form \code{"locus.allele"}.  See
  \code{\link{Genotypes<-,genbinary-method}} for more information.
  }}
}
\details{
  \code{Samples<-} and \code{Loci<-} can only be used to change sample
  and locus names, not to add or remove samples and loci from the
  dataset.

  For slots that require integer values, numerical values used in
  replacement functions will be coerced to integers.  The replacement
  functions also ensure that all slots remain properly indexed.

  The \code{Missing<-} function finds any genotypes with the old missing
  data symbol and changes them to the new missing data symbol, then
  assigns the new symbol to the slot that indicates what the missing
  data symbol is.  \code{Present<-} and \code{Absent<-} work similarly for
  the \code{genbinary} class.

  The \code{Genotype} access and replacement functions deal with
  individual genotypes, which are vectors in the \code{genambig} class.
  The \code{Genotypes} access and replacement functions deal with lists
  of genotypes.

  The \code{PopInfo<-} replacement function also adds elements to
  \code{PopNames(object)} if necessary in order to have names for all of
  the populations.  These will be of the form \code{"Pop"} followed by
  the population number, and can be later edited using
  \code{PopNames<-}.

  The \code{PopNum<-} replacement function first finds all samples in
  the population \code{popname}, and replaces the number in
  \code{PopInfo(object)} for those samples with \code{value}.  It then
  inserts \code{NA} into the original \code{PopNames} slot that
  contained \code{popname}, and inserts \code{popname} into
  \code{PopNames(object)[value]}.  If this results in two populations
  being merged into one, a message is printed to the console.
}
\value{
  \code{PopInfo}, \code{PopNames}, \code{Missing}, \code{Description},
  \code{Usatnts}, \code{Ploidies} and \code{Genotypes} simply return the
  contents of the slots of the same names (or in the case of
  \code{Ploidies}, \code{object@Ploidies@pld} is returned). \code{Samples} and
  \code{Loci} return character vectors taken from the \code{names} of
  other slots (\code{PopInfo} and \code{Usatnts},
  respectively; the initialization and replacement methods ensure that
  these slots are always named according to samples and
  loci). \code{PopNum} returns an integer vector indicating the
  population number(s) of the population(s) named in \code{popname}.
  \code{Genotype} returns a single genotype for a given sample and
  locus, which is a vector whose exact form will depend on the class of
  \code{object}.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{deleteSamples}}, \code{\link{deleteLoci}},
\code{\link{viewGenotypes}}, \code{\link{editGenotypes}},
\code{\link{isMissing}}, \code{\link{estimatePloidy}},
\code{\link[=merge]{merge,gendata,gendata-method}}, \code{\linkS4class{gendata}}
}
\examples{
# create a new genambig (subclass of gendata) object to manipulate
mygen <- new("genambig", samples=c("a", "b", "c"), loci=c("locG",
"locH"))

# retrieve the sample and locus names
Samples(mygen)
Loci(mygen)

# change some of the sample and locus names
Loci(mygen) <- c("lG", "lH")
Samples(mygen)[2] <- "b1"

# describe the dataset
Description(mygen) <- "Example dataset for documentation."

# name some populations and assign samples to them
PopNames(mygen) <- c("PopL", "PopK")
PopInfo(mygen) <- c(1,1,2)
# now we can retrieve samples by population
Samples(mygen, populations="PopL")
# we can also adjust the numbers if we want to make them
# match another dataset
PopNum(mygen, "PopK") <- 3
PopNames(mygen)
PopInfo(mygen)
# change the population identity of just one sample
PopInfo(mygen)["b1"] <- 3

# indicate that both loci are dinucleotide repeats
Usatnts(mygen) <- c(2,2)

# indicate that all samples are tetraploid
Ploidies(mygen) <- 4
# or
Ploidies(mygen) <- rep(4, times = length(Samples(mygen)) * length(Loci(mygen)))
# actually, one sample is triploid
Ploidies(mygen)["c",] <- 3
# view ploidies
Ploidies(mygen)

# view the genotype array as it currently is: filled with missing
# values
Genotypes(mygen)
# fill in the genotypes
Genotypes(mygen, loci="lG") <- list(c(120, 124, 130, 136), c(122, 120),
                                    c(128, 130, 134))
Genotypes(mygen, loci="lH") <- list(c(200, 202, 210), c(206, 208, 210,
                                                        214),
                                    c(208))
# genotypes can also be edited or retrieved by sample
Genotypes(mygen, samples="a")
# fix a single genotype
Genotype(mygen, "a", "lH") <- c(200, 204, 210)
# retrieve a single genotype
Genotype(mygen, "c", "lG")

# change a genotype to being missing
Genotype(mygen, "c", "lH") <- Missing(mygen)
# show the current missing data symbol
Missing(mygen)
# an example of genotypes where one contains the missing data symbol
Genotypes(mygen, samples="c")
# change the missing data symbol
Missing(mygen) <- as.integer(-1)
# now look at the genotypes
Genotypes(mygen, samples="c")

}
\keyword{ methods }

