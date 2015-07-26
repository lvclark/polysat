\name{editGenotypes}
\alias{editGenotypes}

\title{
Edit Genotypes Using the Data Editor
}
\description{
The genotypes from an object of one of the subclasses of \code{gendata}
are converted to a data frame (if necessary), then displayed in the data
editor.  After
the user makes the desired edits and closes the data editor window, the
new genotypes are written to the \code{gendata} object and the object is
returned.
}
\usage{
editGenotypes(object, maxalleles = max(Ploidies(object)),
              samples = Samples(object), loci = Loci(object))
}

\arguments{
  \item{object}{
    An object of the class \code{genambig} or \code{genbinary}.  Contains
    the genotypes to be edited.
}
  \item{maxalleles}{
    Numeric.  The maximum number of alleles found in any given genotype.
    The method
for \code{genambig} requires this information in order to determine how
many columns to put in the data frame.
}
  \item{samples}{
Character or numeric vector indicating which samples to edit.
}
  \item{loci}{
Character or numeric vector indicating which loci to edit.
}
}
\details{
The method for \code{genambig} lists sample and locus names in each row
  in order to identify the genotypes.  However, only the alleles
  themselves should be edited.  NA values and duplicate alleles in the
  data editor will be
  omitted from the genotype vectors that are written back to the
  \code{genambig} object.
}
\value{
An object identical to \code{object} but with edited genotypes.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{viewGenotypes}}, \code{\link{Genotype<-}},
\code{\link{Genotypes<-}}
}
\examples{
if(interactive()){  #this line included for automated checking on CRAN

# set up "genambig" object to edit
mygen <- new("genambig", samples = c("a", "b", "c"),
             loci = c("loc1", "loc2"))
Genotypes(mygen, loci="loc1") <- list(c(133, 139, 142),
                                      c(130, 136, 139, 145),
                                      c(136, 142))
Genotypes(mygen, loci="loc2") <- list(c(202, 204), Missing(mygen),
                                      c(200, 206, 208))
mygen <- reformatPloidies(mygen, output="one")
Ploidies(mygen) <- 4

# open up the data editor
mygen <- editGenotypes(mygen)

# view the results of your edits
viewGenotypes(mygen)

}
}

\keyword{ manip }

