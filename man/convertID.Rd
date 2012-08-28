\name{convertID}
\alias{convertID}
\title{convert gene IDs by id map}
\description{
For same gene, there are multple gene alias. In order to eliminate the possibility of 
missing any connections, convert the gene symbols to unique gene ids is important. This 
function can convert the gene symbols to unique ids and convert it back according a giving
map.
}
\usage{
	convertID(x,IDsMap,ByName=c("from","to"))
}
\arguments{
	\item{x}{a matrix or dataframe contain the columns to be converted.}
	\item{IDsMap}{a character vector of the identifier map}
	\item{ByName}{the column names to be converted}
}
\value{
	a matrix or dataframe with converted gene IDs
}
\examples{
	data("ce.IDsMap")
	bind<-cbind(from="daf-16", to=c("fkh-7", "hlh-13", "mxl-3", "nhr-3", "lfi-1"))
	convertID(toupper(bind), ce.IDsMap, ByName=c("from", "to"))
}
\keyword{convert}