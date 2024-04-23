\name{fpgpnorm}
\alias{fpgpnorm}
\title{Get GP norms from a NODE file for Full-POS fields}
\description{
Extract the GP norm information from a NODE file for Full-POS fields
}
\usage{fpgpnorm(nd,lev,lines.only=FALSE)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
\item{lev}{integer value indicating the level to get norms from. Level 0 means average
 of 2D and 3D fields, -1 means all levels for 3D fields. Any other value should be
 inside range 1:NFLEVG. If a vector is provided, it should be the whole range 1:NFLEVG.}
\item{lines.only}{logical. If true, only GP norm heading lines are returned.}
}
\details{
Arrays have these dimensions: number of steps, number of levels (1 if levels are not asked
 for) and 3. This last dimension contains values for average, minimum and maximum values
 of the field.
}
\value{
If lines.only is false (default), an array of steps, levels and statistical values.
If true, lines of GP norm for Full-POS fields are returned.
}
\seealso{
fpspnorm, spnorm, gpnorm
}
\examples{
nd = readLines(system.file("NODE.001_01",package="mfnode"))
fp = fpgpnorm(nd,0)
str(fp)
}
