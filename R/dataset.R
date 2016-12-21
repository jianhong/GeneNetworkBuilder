#' example datasets for documentation
#' 
#' @description example.data is a data list of example datasets.
#' There is a dataset example.microarrayData, which is the example of gene expression data of a gene-chip result of \emph{C.elegans}.
#' Dataset example.data$ce.bind is a TF binding matrix of ChIP-chip experiment of \emph{C.elegans}.
#' Dataset example.data$cd.exprData is expression data of a gene-chip result  of \emph{C.elegans}.
#' Dataset example.data$hs.bind is a TF binding matrix of ChIP-chip experiment of \emph{H.sapiens}.
#' Dataset example.data$hs.exprData is expression data of a combination of a gene-chip result and a RNA-SEQ result of \emph{H.sapiens}.
#' 
#' @format dataframe
#' @details
#' The dataset example.microarrayData contains columns: ID, logFC, AveExpr, t, P.Value, adj.P.Val, B, genes and symbols. The columns of ID, logFC and symbols are required by GeneNetworkBuilder.
#' The dataset example.data$hs.bind contains columns: ID, symbols, logFC and P.Value.
#' The dataset example.data$hs.exprData contains columns: from and to.
#'
#' @examples
#' data(example.data)
#' names(example.data)
#' head(example.data$example.microarrayData)
#' head(example.data$ce.bind)
#' head(example.data$ce.exprData)
#' head(example.data$hs.bind)
#' head(example.data$hs.exprData)
#' 
"example.data"
#> [1] "example.data"


#' C.elegns gene name to wormbase identifier map
#'
#' map file for converting gene name or sequence name of 
#' \emph{Caenorhabditis elegans} to wormbase identifier
#'
#' @format character vector
#' @details character vecotr with gene name or sequence name as names and
#'           wormbase identifier as values.
#' @source \url{http://www.wormbase.org/}
#' @examples 
#' data(ce.IDsMap)
#' head(ce.IDsMap)
#' 
"ce.IDsMap"
#> [1] "ce.IDsMap"


#' transcript regulatory map of \emph{Caenorhabditis elegans}
#'
#' transcript regulatory map of \emph{Caenorhabditis elegans}
#'
#' @format datafram
#' @details transcript regulatory map of \emph{Caenorhabditis elegans} is 
#'          generated using databases edgedb and microCosm Targets.
#' @source \url{http://edgedb.umassmed.edu}, \url{http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/}
#' @examples 
#' data(ce.interactionmap)
#' head(ce.interactionmap)
#' 
"ce.interactionmap"
#> [1] "ce.interactionmap"


#' map file for converting from wormbase identifier to \emph{Caenorhabditis elegans} gene name
#'
#' map file for converting from wormbase identifier to \emph{Caenorhabditis elegans} gene name
#'
#' @format character vector
#' @details character vecotr with wormbase identifier as names and gene name as values.
#' @source \url{http://www.wormbase.org/}
#' @examples 
#' data(ce.mapIDs)
#' head(ce.mapIDs)
#' 
"ce.mapIDs"
#> [1] "ce.mapIDs"


#' micro RNA of \emph{Caenorhabditis elegans}
#'
#' micro RNA of \emph{Caenorhabditis elegans}
#'
#' @format dataframe
#' @details The first column is wormbase identifier. And the second column is miRNA names.
#' @source \url{http://www.mirbase.org/}
#' @examples 
#' data(ce.miRNA.map)
#' head(ce.miRNA.map)
#' 
"ce.miRNA.map"
#> [1] "ce.miRNA.map"


#' map file for converting gene name or sequence name of \emph{Homo sapiens} to Entrez identifier
#'
#' map file for converting gene name or sequence name of \emph{Homo sapiens} to Entrez identifier
#'
#' @format character vector
#' @details character vecotr with gene name as names and Entrez identifier as values.
#' @examples 
#' data(hs.IDsMap)
#' head(hs.IDsMap)
#' 
"hs.IDsMap"
#> [1] "hs.IDsMap"


#' transcript regulation map of \emph{Homo sapiens}
#'
#' transcript regulation map of \emph{Homo sapiens}
#'
#' @format datafram
#' @details transcript regulatory map of \emph{Homo sapiens} is generated using 
#'          databases FANTOM, mirGen and microCosm Targets.
#' @source \url{http://fantom.gsc.riken.jp/5/}, 
#'         \url{http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/},
#'         \url{http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php}
#' @examples 
#' data(hs.interactionmap)
#' head(hs.interactionmap)
#' 
"hs.interactionmap"
#> [1] "hs.interactionmap"


#' map file for converting from Entrez identifier to \emph{Homo sapiens} gene name
#' 
#' map file for converting from Entrez identifier to \emph{Homo sapiens} gene name
#'
#' @format character vector
#' @details character vecotr with Entrez identifier as names and gene name as values.
#' @examples 
#' data(hs.mapIDs)
#' head(hs.mapIDs)
#' 
"hs.mapIDs"
#> [1] "hs.mapIDs"


#' micro RNA of \emph{Homo sapiens}
#'
#' micro RNA of \emph{Homo sapiens}
#'
#' @format dataframe
#' @details The first column is entrez identifier. And the second column is miRNA names.
#' @source \url{http://www.mirbase.org/}
#' @examples 
#' data(hs.miRNA.map)
#' head(hs.miRNA.map)
#' 
"hs.miRNA.map"
#> [1] "hs.miRNA.map"
