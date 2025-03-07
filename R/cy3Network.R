#' visualize the network by Cy3
#' @description
#' Using RCy3 to visualize the network.
#' 
#' @param gR an object of \link[graph:graphNEL-class]{graphNEL}
#' @param ... parameters will be passed to \link[RCy3]{createNetworkFromGraph}
#' @param stringify Run STRINGify or not
#' @param species if stringify is TRUE, teh species will be used to retreive
#' the nodes and edges properties.
#' @param style The default style when create the network
#' @param widths The link width range.
#' @return The network SUID.
#' @importFrom RCy3 createNetworkFromGraph getVisualStyleNames setVisualStyle
#' setNodeSizeMapping setNodeBorderColorBypass setNodeColorBypass 
#' setEdgeLineWidthMapping commandsRun getNetworkSuid getTableColumns
#' lockNodeDimensions getInstalledApps
#' @export
#' @examples
#' data("ce.miRNA.map")
#' data("example.data")
#' data("ce.interactionmap")
#' data("ce.IDsMap")
#' sifNetwork<-buildNetwork(example.data$ce.bind, ce.interactionmap, level=2)
#' cifNetwork<-filterNetwork(rootgene=ce.IDsMap["DAF-16"], sifNetwork=sifNetwork, 
#'   exprsData=uniqueExprsData(example.data$ce.exprData, "Max", condenseName='logFC'),
#'   mergeBy="symbols",
#'   miRNAlist=as.character(ce.miRNA.map[ , 1]), tolerance=1)
#' gR<-polishNetwork(cifNetwork)
#' if(interactive()){
#'   cy3Network(gR)
#' }
cy3Network <- function(gR = graphNEL(), ...,
                       stringify=FALSE, species = 'Homo sapiens',
                       style='Marquee', widths = c(0.25, 5)){
  stopifnot(is(gR, 'graphNEL'))
  dots <- list(...)
  if('base.url' %in% names(dots)){
    base.url <- dots[["base.url"]]
  }else{
    base.url <- 'http://127.0.0.1:1234/v1'
  }
  availableStyles <- getVisualStyleNames(base.url=base.url)
  style <- match.arg(style, choices = availableStyles)
  network <- createNetworkFromGraph(gR, ...)
  if(stringify){
    installedapps <- getInstalledApps(base.url = base.url)
    if(!grepl('stringApp', installedapps)){
      warning('Please install STRING app first!')
    }
    ## doc https://github.com/RBVI/stringApp/blob/master/doc_automation.md#compound-query
    string.cmd = paste0('string stringify column="name" networkNoGui="current" species="', species, '"')
    commandsRun(string.cmd, base.url = base.url)
    
    network <- getNetworkSuid(base.url = base.url)
    stringstyle <- "STRING - From graph"
    if(!stringstyle %in% getVisualStyleNames(base.url = base.url)){
      stop('Default string style is not in list.',
           'Please select correct one from the output of getVisualStyleNames()')
    }else{
      style <- stringstyle
    }
  }else{
    setVisualStyle(style, network = network, base.url = base.url)
    lockNodeDimensions(TRUE, style.name = style, base.url = base.url)
    setNodeSizeMapping('size', style.name = style, network = network, base.url = base.url)
  }
  nodeData <- getTableColumns(table = 'node', network = network, base.url = base.url)
  borderColor <- graph::nodeRenderInfo(gR, 'col')
  if(length(borderColor)>0){
    if(stringify){
      N <- nodeData$name[match(names(borderColor), nodeData$`query term`)]
    }else{
      N <- names(borderColor)
    }
    borderColor[is.na(borderColor)] <- '#000000'
    setNodeBorderColorBypass(N[!is.na(N)], 
                             borderColor[!is.na(N)],
                             network = network,
                             base.url = base.url)
  }
  color <- graph::nodeRenderInfo(gR, 'fill')
  if(length(color)>0){
    if(stringify){
      N <- nodeData$name[match(names(color), nodeData$`query term`)]
    }else{
      N <- names(color)
    }
    color[is.na(color)] <- '#FFFFFF'
    setNodeColorBypass(N[!is.na(N)],
                       new.colors=color[!is.na(N)],
                       network = network,
                       base.url = base.url)
  }
  setEdgeLineWidthMapping('weight',
                          style.name = style,
                          widths = widths,
                          network = network,
                          base.url = base.url)
  return(network)
}
