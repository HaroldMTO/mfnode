\name{gpnorm}
\alias{gpnorm}
\title{Array of GP norms}
\description{
Extract to an array the GP norm information from a NODE file (via its lines)
}
\usage{gpnorm(nd,lev=0,tag="NORMS AT (START|NSTEP|END) CNT4",gpin="\\w+.*",
 gpout=character,lines.only=FALSE,abbrev=TRUE)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
\item{lev}{integer value giving the level to get GP norm information from.
 If lev=0, the average of levels is taken. If lev > 0, level lev is taken.
 Alternatively, lev can be the vector of NFLEVG levels (ie all levels or 1:NFLEVG).
}
\item{tag}{regular expression matching the heading line for GP norm lines}
\item{gpin}{regular expression matching field names to include in GP norm lines}
\item{gpout}{regular expression matching field names to exclude from GP norm lines}
\item{lines.only}{logical. If true, only GP norm heading lines are returned.}
\item{abbrev}{logical. If true, shorten dimnames 3, that is, names of variables.}
}
\details{
Default values are made for classical GFL norm lines printed from CNT4.
Array has these dimensions: number of steps, number of levels (=1 if 3D average), 3
 and number of GP fields.
The third dimension, 3, contains values for average, minimum and maximum values of the
 field.
}
\value{
If lines.only is false (default), an array of GP norms.
If true, heading lines of GP norms are returned.
}
\seealso{
spnorm, fpgpnorm, fpspnorm
}
\examples{
nd = readLines(system.file("data/NODE.001_01",package="mfnode"))
gfl0 = gpnorm(nd,lev=0)
str(gfl0)
gmv0 = gpnorm(nd,lev=0,tag="gpnorm gmvt0")
gmv1 = gpnorm(nd,lev=0,tag="gpnorm gmvt1 lag")
str(gmv0)
str(gmv1)
}
