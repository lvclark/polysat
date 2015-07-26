\name{reformatPloidies}
\alias{reformatPloidies}
\title{
Convert Ploidy Format of a Dataset
}
\description{
This function changes the class of the object in the \code{Ploidies}
slot of a \code{"gendata"} object.  (See the four subclasses described
in \code{"\linkS4class{ploidysuper}"}.)  Existing ploidy data can either
be erased or, if possible, used in the new format.
}
\usage{
reformatPloidies(object, output = "collapse", na.rm = FALSE, erase = FALSE)
}
\arguments{
  \item{object}{
A \code{"gendata"} object.
}
  \item{output}{
A character string indicating the desired result of the conversion:
\code{"matrix"} if ploidies should be indexed by both sample and locus,
\code{"sample"} if ploidies should be indexed only by sample,
\code{"locus"} if ploidies should be indexed only by locus, \code{"one"}
if there should be one ploidy for the entire dataset, or
\code{"collapse"} if ploidies should be converted to the simplest
possible format.
}
  \item{na.rm}{
Boolean.  If \code{FALSE}, \code{NA} is treated as a unique ploidy.  If
\code{TRUE}, \code{NA} values are ignored assuming that each sample
and/or locus has only one ploidy otherwise.  This argument is passed
directly to \code{\link{plCollapse}}.
}
  \item{erase}{
Boolean.  If \code{TRUE}, the new \code{Ploidies} slot is simply filled
with \code{NA} instead of existing ploidy values from \code{object}.
}
}
\details{
  This is a versatile function that can accomplish several tasks
  relating to the format of ploidies in the dataset:

  If you wish to change how ploidy is indexed, but don't care about
  keeping any data in the \code{Ploidies} slot, set \code{erase=TRUE}
  and \code{output} to \code{"matrix"}, \code{"sample"}, \code{"locus"},
  or \code{"one"}.

  If you wish to keep ploidy data while moving from a simpler format to
  a more complex format (\emph{i.e.} from \code{"one"} to any other
  format, or from \code{"sample"} or \code{"locus"} to \code{"matrix"}),
  leave \code{erase=FALSE} and set \code{output} to the desired format.
  Existing data will be duplicated to fill out the new format.  For
  example, if ploidies were indexed by sample and you change to matrix
  format, the ploidy that had previously been recorded for each sample
  will be duplicated for each locus.

  If you wish to keep ploidy data while performing any other format
  conversions (\emph{e.g.} \code{"matrix"} to \code{"sample"} or
  \code{"sample"} to \code{"locus"}), the function will check that there
  is one unique ploidy for each sample, locus, or the entire dataset (as
  appropriate), and will produce an error if the conversion cannot be
  done without a loss of information.

  If you wish to keep ploidy data and convert to the simplest possible
  format, set \code{output="collapse"}.  The function will automatically
  determine the simplest format for conversion without loss of data.
  (The \code{read} functions in \pkg{polysat} that take ploidy data from
  the input file use this option.)
}
\value{
A \code{"gendata"} object that is a copy of \code{object} but with the
\code{Ploidies} slot converted to a new class.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{plCollapse}}, \code{"\linkS4class{ploidysuper}"}
}
\examples{
# Make a new "genambig" object for this example
testdata <- new("genambig")
Ploidies(testdata)

## If you need to reformat before you have entered any ploidy
## information:

# convert from matrix to sample format
testdata <- reformatPloidies(testdata, output="sample")
Ploidies(testdata)

## If you have entered ploidy information but realized you can use a
## simpler format:

# Enter some ploidies
Ploidies(testdata)[1] <- 2
Ploidies(testdata)

# Convert from "sample" to "one" with na.rm=TRUE
testdata <- reformatPloidies(testdata, na.rm=TRUE, output="one")
Ploidies(testdata)

## If you change your mind and want to go back to a more complex format

# Convert from "one" to "locus"
testdata <- reformatPloidies(testdata, output="locus")
Ploidies(testdata)
}
\keyword{ manip }

