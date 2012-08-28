\name{buildNetwork}
\alias{buildNetwork}
\title{construct the regulatory network}
\description{
get all the connections of interesting genes from regulatory map.
}
\usage{
buildNetwork(TFbindingTable, interactionmap, level=3)
}
\arguments{
\item{TFbindingTable}{a matrix or data.frame with interesting genes. Column names must be 'from', 'to'}
\item{interactionmap}{Transcription regulatory map. Column names of interactionmap must
 be 'from','to'}
\item{level}{Depth of node path}
}
\value{
 a dataframe or matrix of all the connections of interesting genes
}
\examples{
data("ce.interactionmap")
data("example.data")
xx<-buildNetwork(example.data$ce.bind, ce.interactionmap, level=2)
}
\keyword{network}