\name{generateCosmic}
\alias{generateCosmic}
\title{
Generate GRanges COSMIC objects
}
\description{
Generate GRanges COSMIC objects from a data.frame read in after
download.  The end user is not meant to use this function.  It is here
to document the approach used to create the data.
}
\usage{
generateCosmic(dat)
}
\arguments{
  \item{dat}{A data.frame after reading a cosmic Mutant Export file
}
}
\value{
\item{b37 }{A GRanges object for build 37}
\item{b36 }{A GRanges object for build 36}

}
\author{
Sean Davis <seandavi@gmail.com>
}
\seealso{
See Also as \code{\link{cosmicGRCh37}}, \code{\link{cosmicNCBI36}}
}
\examples{
generateCosmic
}
\keyword{ IO }
\keyword{ datasets }
