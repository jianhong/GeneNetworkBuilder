\name{uniqueExprsData}
\alias{uniqueExprsData}
\title{unique the microarray data}
\description{
get unique the microarray data for each gene id.
}
\usage{
uniqueExprsData(exprsData, method, condenseName='logFC')
}
\arguments{
\item{exprsData}{dataset of expression comparison data}
\item{method}{method must be Max, Median or Min}
\item{condenseName}{column names to be condensed}
}
\value{
 a dataframe of expression data without duplicates
}
\examples{
data("example.data")
example.microarrayData<-uniqueExprsData(example.data$example.microarrayData,
						"Max", condenseName='logFC')
}
\keyword{network}