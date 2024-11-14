\name{scale10}
\alias{scale10}
\title{Data scale}
\description{
Compute the power of 10 that corresponds most to x
}
\usage{scale10(x,scientific=TRUE)}

\arguments{
\item{x}{numeric vector (or matrix, array)}
\item{scientific}{logical indicating whether scale should be a power of 10 (FALSE) 
 or 1000 (TRUE)}
}
\details{
NaN values are removed from input before scaling.
If resulting x data is of length 0, NaN is returned with a warning.
Scale is computed for the absolute maximum (ie. 'max(abs(x),na.rm=TRUE)').
Scale of 0 is 1.
}
\value{
If scientific is TRUE (default), scale of x is given as a power of 1000, 10 if FALSE.
}
\seealso{
intlines
}
\examples{
scale10(1)
scale10(-1,FALSE)
scale10(Inf)
scale10(c(0,NA))
}
