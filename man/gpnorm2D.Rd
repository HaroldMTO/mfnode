\name{gpnorm}
\alias{gpnorm}
\title{List of 2D GP norms}
\description{
Extract to a list the GP norm information for surface fields from a NODE file
}
\usage{gpnorm2D(nd)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
}
\details{
List contains one element per group of surface fields. The elements found for groups
 are arrays of dimension 2. First dimension is for fields, 2nd one is for statistical
 values (ave/min/max).
}
\value{
A list of arrays of GP norms for surface fields in groups.
}
\seealso{
gpnorm2D, spnorm, fpgpnorm, fpspnorm
}
\examples{
nd = readLines(system.file("data/NODE.001_01",package="mfnode"))
surf = gpnorm2D(nd)
str(surf)
}
