#' Save network in various formats
#' @description Save graph into HTML, json or xgmml format.
#' @importFrom rjson toJSON
#' @import htmlwidgets
#' @export
#' @param network output of \link{browseNetwork}
#' @param file Name of the file to save to.
#' @param format type in which graph shall be saved. 
#'               Could be one of HTML, json or XGMML.
#' @param ... Parameter could be used by \link[htmlwidgets]{saveWidget} for HTML
#'            or writeLines for json or \link[XML]{saveXML} for XGMML.
#' @examples
#' data("ce.miRNA.map")
#' data("example.data")
#' data("ce.interactionmap")
#' data("ce.IDsMap")
#' sifNetwork<-buildNetwork(example.data$ce.bind, ce.interactionmap, level=2)
#' cifNetwork<-filterNetwork(rootgene=ce.IDsMap["DAF-16"], sifNetwork=sifNetwork, 
#'                         exprsData=uniqueExprsData(example.data$ce.exprData, "Max", condenseName='logFC'),
#'                         mergeBy="symbols",
#'                         miRNAlist=as.character(ce.miRNA.map[ , 1]), tolerance=1)
#' gR<-polishNetwork(cifNetwork)
#' network <- browseNetwork(gR)
#' exportNetwork(network, "sample.html")
#' @keywords IO
exportNetwork <- function(network, file, 
                          format=c("HTML", "json", "XGMML"), ...){
  format <- match.arg(format)
  stopifnot(inherits(network, c("browseNetwork")))
  if(length(network$x$elements$nodes)<1){
    stop("There is no node in the network.")
  }
  return(invisible(switch(format,
                          "HTML"={saveWidget(network, file, ...)},
                          "json"={writeLines(toJSON(network$x), con = file, ...)},
                          "XGMML"={saveXGMML(network$x, file = file, ...)})))
}

#' Save network as xgmml
#' @description Save graph into xgmml format.
#' @param network output of \link{browseNetwork}
#' @param file Name of the file to save to.
#' @param ... Parameter could be used by saveXML
#' @importFrom XML xmlNode append.xmlNode saveXML
#' @keywords IO
#' 
saveXGMML <- function(network, file, ...){
  top <- xmlNode("graph", 
                 attrs=c("label"="GeneNetworkBuilder",
                         "xmlns:dc"="http://purl.org/dc/elements/1.1/", 
                         "xmlns:xlink"="http://www.w3.org/1999/xlink", 
                         "xmlns:rdf"="http://www.w3.org/1999/02/22-rdf-syntax-ns#", 
                         "xmlns:cy"="http://www.cytoscape.org", 
                         "xmlns"="http://www.cs.rpi.edu/XGMML",
                         "directed"="1"))
  attr <- names(network$elements$nodes[[1]]$data)
  attr.cl <- sapply(network$elements$nodes[[1]]$data, class)
  attr.map <- c("character"="string", "numeric"="real", 
                "logical"="boolean", "integer"="interger")
  nodes <- lapply(network$elements$nodes, function(.ele){
    n <- xmlNode("node",
                 attrs=c("label"=.ele$data$label,
                         "id"=.ele$data$id))
    for(i in c("shared name", "name")){
      n <- append.xmlNode(n,
                          xmlNode("att",
                                  attrs=c("name"=i,
                                          "value"=.ele$data$label,
                                          "type"="string")))
    }
    for(att in attr){
      n <- append.xmlNode(n,
                          xmlNode("att", 
                                  attrs=c("name"=att,
                                          "type"=attr.map[attr.cl[att]],
                                          "value"=.ele$data[[att]])))
    }
    n <- append.xmlNode(n,
                        xmlNode("graphics",
                                attrs=c("outline"=.ele$data$borderColor,
                                        "width"=2,
                                        "x"=.ele$data$nodeX,
                                        "h"=.ele$data$size,
                                        "y"=.ele$data$nodeY,
                                        "w"=.ele$data$size,
                                        "type"="ELLIPSE",
                                        "fill"=.ele$data$fill)))
    n
  })
  for(n in nodes){
    top <- append.xmlNode(top, n)
  }
  edges <- lapply(network$elements$edges, function(.ele){
    e <- xmlNode("edge",
                 attrs=c("label"=.ele$data$label,
                         "source"=.ele$data$source,
                         "target"=.ele$data$target))
    e <- append.xmlNode(e, 
                        xmlNode("att",
                                attrs=c("name"="weight",
                                        "type"="real",
                                        "value"=.ele$data$weight)))
    e
  })
  for(e in edges){
    top <- append.xmlNode(top, e)
  }
  saveXML(top, file=file, ...)
}