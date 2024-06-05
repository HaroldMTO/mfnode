\name{spnorm}
\alias{spnorm}
\title{Parse lines to integer}
\description{
Parse lines of integer values into an integer vector
}
\usage{intlines(nd)}

\arguments{
\item{nd}{character vector of lines read from a NODE file}
}
\details{
Any integer value in the character lines passed as nd are parsed and returned.
Integer values can appear after some non numeric part in the lines.
For better parsing, character lines can be processed for some subtitution.
}
\value{
An integer vector of the values parsed.
}
\seealso{
intlines
}
\examples{
nd = readLines(system.file("data/NODE.001_01",package="mfnode"))
ind = grep("SETA=",nd,ignore.case=TRUE)
x = intlines(nd[ind])
ind = grep("ARRAY.+ALLOCATED",nd)
y = intlines(sub("ARRAY.+ALLOCATED","",nd[ind]))
str(x)
str(y)
}
