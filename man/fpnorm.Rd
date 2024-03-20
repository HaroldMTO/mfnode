\name{fpnorm}
\alias{fpnorm}
\title{List of arrays of GP norms for Full-POS fields}
\description{
Extract to a list of arrays the GP norm information for Full-POS fields from a NODE file
 (via its lines)
}
\usage{fpnorm(nd,lines.only=FALSE)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
\item{lines.only}{logical. If true, only GP norm heading lines are returned.}
}
\details{
2D and 3D fields appear in list as one element for each field. 3D fields have one
dimension more as 2D fields.
Arrays have these dimensions: number of steps, number of levels (only if 3D field)
 and 3.
The last dimension, 3, contains values for average, minimum and maximum values of the
 field.
}
\value{
If lines.only is false (default), a list composed of one array for each field.
If true, lines of GP norm for Full-POS fields are returned.
}
\seealso{
spnorm, gpnorm
}
\examples{
nd = readLines(system.file("NODE.001_01",package="mfnode"))
fpl = fpnorm(nd)
str(fpl)
}
