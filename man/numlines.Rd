\name{spnorm}
\alias{spnorm}
\title{Parse lines to numeric}
\description{
Parse lines of numeric values into a numeric vector
}
\usage{numlines(nd)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
}
\details{
Numeric values can appear after some non numeric part in the lines.
}
\value{
A numeric vector of the values parsed.
}
\seealso{
intlines
}
\examples{
nd = readLines(system.file("data/NODE.001_01",package="mfnode"))
ind = grep("A and B",nd,ignore.case=TRUE)
nflevg = getvar("NFLEVG",nd)
x = numlines(nd[ind+1+seq(nflevg)])
str(x)
}
