\name{spnorm}
\alias{spnorm}
\title{Array of SP norms}
\description{
Extract to an array the SP norm information from a NODE file
}
\usage{spnorm(nd,lev=0,tag="NORMS AT (ART|NSTEP|END) CNT4",abbrev=TRUE)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
\item{lev}{integer value (or vector) giving the level to get GP norm information from.
 If lev=0, the average of levels is got. If lev > 0, level lev is got. If lev is a vector,
 it must be the range of NFLEVG levels (ie all levels).}
\item{tag}{regular expression matching the heading line for SP norm lines}
\item{abbrev}{logical. If true, shorten dimnames 3, that is, names of variables.}
}
\details{
Function gpnormhead prepares a call to gpnorm. Default values are made for GFL norm lines
 printed from CNT4.
Array has these dimensions: number of steps, number of levels (=1 if 3D average),
 number of fields.
}
\value{
An array of SP norms with steps, levels and variables found.
}
\seealso{
gpnorm, gpnorm2D, fpspnorm, fpgpnorm
}
\examples{
nd = readLines(system.file("data/NODE.001_01",package="mfnode"))
sp0 = spnorm(nd,lev=0)
str(sp0)
sp1 = spnorm(nd,lev=0,tag="spnorm t1si")
str(sp1)
}
