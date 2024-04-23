\name{gpnorm}
\alias{gpnorm}
\title{List of 2D GP norms}
\description{
Extract to a list the GP norm information for surface fields from a NODE file
}
\usage{gpnorm2D(nd,lines.only=FALSE)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
\item{lines.only}{logical. If true, only GP norm heading lines are returned.}
}
\details{
List contains one element per group of surface fields. The elements found for groups
 are arrays of dimension 2. First dimension is for fields, 2nd one is for statistical
 values (ave/min/max).
}
\value{
If lines.only is false (default), a list of arrays of GP norms for surface fields in
 groups.
If true, heading lines of groups of GP norms are returned.
}
\seealso{
gpnorm2D, spnorm, fpgpnorm, fpspnorm
}
\examples{
nd = readLines(system.file("data/NODE.001_01",package="mfnode"))
surf = gpnorm2D(nd)
str(surf)
}
