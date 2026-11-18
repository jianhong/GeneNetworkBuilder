#' unique the microarray data
#' 
#' @description get unique the microarray data for each gene id.
#' 
#' @param exprsData dataset of expression comparison data
#' @param method method must be Max, Median or Min
#' @param condenseName column names to be condensed
#' 
#' @return a dataframe of expression data without duplicates
#' @keywords network
#' @examples 
#' data("example.data")
#' example.microarrayData<-uniqueExprsData(example.data$example.microarrayData,
#'                                         method="Max", condenseName='logFC')
#'
#' @export
#' @importFrom plyr . ddply
#' @importFrom graphics symbols
#' @importFrom stats median
#' 
uniqueExprsData<-function(exprsData, method='Max', condenseName='logFC'){
    if(!(method %in% c("Max", "Median", "Min"))){
        stop("method must be Max, Median or Min")
    }
    if(!checkCName("symbols", exprsData)){
        stop("symbols is not a valide colname of exprsData")
    }
    if(!checkCName(condenseName, exprsData)){
        stop(paste(condenseName," is not a valide colname of exprsData"))
    }
    if(!is.data.frame(exprsData)){
        exprsData<-as.data.frame(exprsData)
    }
    if(!is.numeric(exprsData[ , condenseName])){
        stop(paste("class of", condenseName, "is not a numeric column"))
    }
    exprsData<-switch(method,
                      Max   =plyr::ddply(exprsData, plyr::.(symbols), 
                                         getMax, 
                                         condenseName),
                      Median=plyr::ddply(exprsData, plyr::.(symbols), 
                                         getMedian, 
                                         condenseName),
                      Min   =plyr::ddply(exprsData, plyr::.(symbols),
                                         getMin, 
                                         condenseName)
                           )
    exprsData
}

#' convert gene IDs by id map
#' @description For same gene, there are multple gene alias. 
#' In order to eliminate the possibility of missing any connections, 
#' convert the gene symbols to unique gene ids is important. 
#' This function can convert the gene symbols to unique ids and 
#' convert it back according a giving map.
#' @param x a matrix or dataframe contain the columns to be converted.
#' @param IDsMap a character vector of the identifier map
#' @param ByName the column names to be converted
#' @return a matrix or dataframe with converted gene IDs
#' @examples 
#' data("ce.IDsMap")
#' bind<-cbind(from="daf-16", to=c("fkh-7", "hlh-13", "mxl-3", "nhr-3", "lfi-1"))
#' convertID(toupper(bind), ce.IDsMap, ByName=c("from", "to"))
#' @keywords convert
#' @export
#' 
convertID<-function(x, IDsMap, ByName=c("from", "to")){
    if((!is.character(IDsMap)) | (is.null(IDsMap))){
        stop("invalide IDsMap")
    }
    for(i in 1:length(ByName)){
        if(!checkCName(ByName[i],x)){
            stop(paste(ByName[i],"is not a valide colname of x"))
        }
        x[,ByName[i]]<-IDsMap[as.character(x[,ByName[i]])]
    }
    x
}

#' construct the regulatory network
#' @description Get all the connections of interesting genes from regulatory map.
#' @param TFbindingTable a matrix or data.frame with interesting genes. 
#'                       Column names must be 'from', 'to'
#' @param interactionmap Transcription regulatory map. 
#'                       Column names of interactionmap must be 'from','to'
#' @param level Depth of node path
#' 
#' @return a dataframe or matrix of all the connections of interesting genes
#' @keywords network
#' @examples 
#' data("ce.interactionmap")
#' data("example.data")
#' xx<-buildNetwork(example.data$ce.bind, ce.interactionmap, level=2)
#' @export

buildNetwork<-function(TFbindingTable, interactionmap, level=3){
    checkMap(interactionmap, TFbindingTable)
    if(level>0){
        y<-interactionmap[interactionmap[ , "from"] %in% unique(as.character(TFbindingTable[ , "to"])), ,drop=FALSE]
        y<-unique(y)
        z<-y[!(y[,"to"] %in% TFbindingTable[,"to"]), , drop=FALSE]
        nrow1<-nrow(TFbindingTable)
        TFbindingTable<-rbind(TFbindingTable, y)
        TFbindingTable<-unique(TFbindingTable)
        level<-level-1
        if(level>0){
            nrow2<-nrow(TFbindingTable)
            if(nrow2>nrow1){
                y<-buildNetwork(z, interactionmap, level)
                TFbindingTable<-rbind(TFbindingTable, y)
            }
            TFbindingTable<-unique(TFbindingTable)
        }
    }
    TFbindingTable
}

#' filter the regulatory network table by expression profile
#' @description verify every nodes in the regulatory network by expression profile
#' @param rootgene name of root gene. It must be the ID used in xx regulatory network
#' @param sifNetwork Transcription regulatory network table. 
#'                   Column names of xx must be 'from','to'
#' @param exprsData dataset of expression comparison data, 
#'                  which should contain column logFC and column given by exprsDataByName
#' @param mergeBy The column name contains ID information used to merge with 
#'                'to' column of sifNetwork in exprsData
#' @param miRNAlist vector of microRNA ids.
#' @param remove_miRNA remove miRNA from the network or not. 
#'                     Bool value, TRUE or FALSE
#' @param tolerance maximum number of unverified nodes in each path
#' @param cutoffPVal cutoff p value of valid differential expressed gene/miRNA
#' @param cutoffLFC cutoff log fold change value of a valid differential 
#'                  expressed gene/miRNA
#' @param minify Only keep the best path if multiple paths exists for single node? 
#'               Bool value, TRUE or FALSE
#' @param miRNAtol take miRNA expression into account for tolerance calculation. 
#'                 Bool value, TRUE or FALSE
#' @return a dataframe of filtered regulatory network by expression profile
#' @import Rcpp
#' @useDynLib GeneNetworkBuilder
#' @export
#' 
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
#' @keywords network
#' 
filterNetwork<-function(rootgene, sifNetwork, exprsData, mergeBy="symbols", miRNAlist, remove_miRNA=FALSE,
                    tolerance=0, cutoffPVal=0.01, cutoffLFC=0.5, minify=TRUE, miRNAtol=FALSE)
{
    checkMCName(sifNetwork)
    if(!missing(miRNAlist)){
        if(!is.vector(miRNAlist)){
            stop("miRNAlist should be a vector")
        }
    }
    if(!checkCName(mergeBy, exprsData)){
        stop(paste(mergeBy, "is not a column name of exprsData"))
    }
    if(!checkCName("logFC", exprsData)){
        stop("logFC is not a column name of exprsData")
    }
    if(!is.numeric(exprsData[ , "logFC"])){
        stop("class of exprsData[ , \"logFC\"] is not a numeric column")
    }
    if(!checkCName("P.Value", exprsData)){
        stop("P.Value is not a column name of exprsData")
    }
    if(!is.numeric(exprsData[ , "P.Value"])){
        stop("class of exprsData[ , \"P.Value\"] is not a numeric column")
    }
    if(!is.numeric(cutoffLFC)){
        stop("cutoffLFC is not a numeric")
    }
    if(!is.numeric(cutoffPVal)){
        stop("cutoffPVal is not a numeric")
    }
    if(any(duplicated(exprsData[,mergeBy]))){
        stop("expresData has multiple logFC for same ID. Please try ?uniqueExprsData")
    }
    if(!is.logical(minify)){
        stop("minify is not a logical")
    }
    if(!is.logical(miRNAtol)){
        stop("miRNAtol is not a logical")
    }
    tolerance<-round(tolerance)
    cifNetwork<-merge(sifNetwork, exprsData,
                      by.x="to", by.y=mergeBy, all.x=TRUE)
    
    if(missing(rootgene)){## unrooted
      rootgene <- NA
    }
    if(is.na(rootgene) | is.null(rootgene)){## unrooted
      ## create a fake rootgene
      ## the fake rootgene will connect to all the genes with filtered conditions
      rootgene <- "NANANA"
      fakeroot <- TRUE
      if(rootgene %in% c(cifNetwork$from, cifNetwork$to)){
        rootgene <- make.names(c(cifNetwork$from, cifNetwork$to, rootgene))
        rootgene <- rootgene[length(rootgene)]
      }
      rootlogFC <- cutoffLFC + 1
      cifNetwork_filtered <- cifNetwork[!is.na(cifNetwork[, "logFC"]) &
                                          !is.na(cifNetwork[, "P.Value"]),
                                        , drop=FALSE]
      cifNetwork_filtered <- cifNetwork_filtered[
        abs(cifNetwork_filtered[, "logFC"])>=cutoffLFC &
        cifNetwork_filtered[, "P.Value"]<=cutoffPVal,
        , drop=FALSE]
      cifNetwork_filtered <- 
        cifNetwork[cifNetwork$to %in% cifNetwork_filtered$from, , drop=FALSE]
      cifNetwork_filtered$from <- rootgene
      cifNetwork_filtered <- unique(cifNetwork_filtered)
      cifNetwork <- rbind(cifNetwork, cifNetwork_filtered)
    }else{
      fakeroot <- FALSE
      rootlogFC<-exprsData[exprsData[ , mergeBy] == rootgene, "logFC"]
      rootlogFC<-rootlogFC[!is.na(rootlogFC)]
      rootlogFC<-ifelse(length(rootlogFC) < 1, 0.0, rootlogFC[1])
    }
##   convert NA to 0 for logFC
    cifNetwork.logFC<-cifNetwork[,"logFC"]
    cifNetwork.logFC[is.na(cifNetwork.logFC)]<-0.0
    cifNetwork.pValue<-cifNetwork[,"P.Value"]
    cifNetwork.pValue[is.na(cifNetwork.pValue)]<-1.0
##   label microRNA
    cifNetwork$miRNA<-FALSE
    cifNetwork$dir<-2
    if(!missing(miRNAlist)){
        if(length(miRNAlist)>0){
            cifNetwork$miRNA<-ifelse(cifNetwork$to %in% miRNAlist, TRUE, FALSE)
            cifNetwork$dir<-ifelse(cifNetwork$from %in% miRNAlist, 0, 2)
        }
    }
##   remove micorRNA
    if(remove_miRNA){
        cifNetwork.logFC<-cifNetwork.logFC[!cifNetwork$miRNA]
        cifNetwork.pValue <- cifNetwork.pValue[!cifNetwork$miRNA]
        cifNetwork<-cifNetwork[!cifNetwork$miRNA, ]
    }
    
    cifNetwork.list <- .Call("filterNodes",
                         as.character(cifNetwork$from), 
                         as.character(cifNetwork$to), 
                         cifNetwork$miRNA, 
                         cifNetwork.logFC,
                         cifNetwork.pValue,
                         cifNetwork$dir,
                         nrow(cifNetwork), 
                         rootgene[1],
                         rootlogFC[1],
                         tolerance[1],
                         minify[1],
                         miRNAtol[1],
                         cutoffLFC[1],
                         cutoffPVal[1]
                         )
    cifNetwork.list <- do.call(rbind, lapply(names(cifNetwork.list),
                                           function(.name, .ele){
                                              if(length(.ele[[.name]])>0){
                                                cbind(from=.ele[[.name]], to=.name)
                                              }else{
                                                cbind(from=NA, to=.name)
                                              }
                                            },
                                           cifNetwork.list)
                            )
    if(minify){
      cifNetwork <- merge(cifNetwork, cifNetwork.list)
    }else{
      cifNetwork_s <- merge(cifNetwork, cifNetwork.list)
      nodes <- unique(c(cifNetwork_s$from, cifNetwork_s$to))
      cifNetwork <- cifNetwork[cifNetwork$from %in% nodes &
                                 cifNetwork$to %in% nodes, , drop=FALSE]
    }
    if(fakeroot){
      cifNetwork$from[cifNetwork$from==rootgene] <- NA
    }
    cifNetwork
}

#' generate an object of grahpNEL to represent the regulation network
#' @description generate an object of grahpNEL to represent the regulation network. 
#' Each node will has three attributes: size, borderColor and fill. 
#' The size will be mapped to the length of its edges. The node fill color will 
#' be mapped to logFC.
#' @param cifNetwork dataframe used to draw network graph. column names of 
#'                   cifNetwork must contain 'from', 'to', 'logFC' and 'miRNA'
#' @param nodeData The node data. If it is not provide, node data will be 
#' retrieved from cifNetwork for the 'to' nodes.
#' @param nodesDefaultSize nodes default size
#' @param nodecolor a character vector of color set. 
#'                  The node color will be mapped to color set by log fold change.
#'                  Or the column names for the colors.
#' @param nodeBg background of node
#' @param nodeBorderColor a list of broder node color set. 
#'                        nodeBorderColor's element must be gene and miRNA
#' @param edgeWeight the weight of edge. It can be a column name of cifNetwork.
#' @param edgelwd the default width of edge. If edgeWeight is set, the edgelwd 
#' will be mapped to the edgeWeight.
#' @param ... any parameters can be passed to \link[graph:settings]{graph.par}
#' @return An object of graphNEL class of the network
#' @import graph
#' @importFrom grDevices colorRampPalette
#' @importFrom methods new
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
#' ##	browseNetwork(gR)
#' @keywords network
#' 
polishNetwork<-function(cifNetwork,
                        nodeData,
                        nodesDefaultSize=48,
                        nodecolor=colorRampPalette(c("green", "yellow", "red"))(5),
                        nodeBg="white",
                        nodeBorderColor=list(gene='darkgreen',miRNA='darkblue'), 
                        edgeWeight=NA, edgelwd=0.5, ...)
{
  cname<-c("from", "to")
  if(!is.data.frame(cifNetwork)){
    stop("cifNetwork should be a data.frame")
  }
  
  cifNetwork<-cifNetwork[!duplicated(cifNetwork[,cname]), ]
  edge<-cifNetwork[cifNetwork$from!="" & cifNetwork$to!="", cname]
  node<-c(as.character(unlist(edge)))
  node<-node[!is.na(node)]
  node<-node[node!=""]
  node<-unique(node)
  if(length(node) <= 1){
    stop("Can not built network for the inputs. Too less connections.")
  }
  
  if(missing(nodeData)){
    if(length(intersect(c("from", "to", "logFC"), colnames(cifNetwork)))<3){
      stop("colnames of cifNetwork must contain 'from', 'to', and 'logFC'");
    }
    nodeData <- cifNetwork[, !colnames(cifNetwork) %in% c("from", "dir"),
                           drop=FALSE]
    nodeData <- nodeData[!duplicated(nodeData[, "to"]), , drop=FALSE]
    rownames(nodeData) <- nodeData[, 'to', drop=TRUE]
    nodeData$to <- NULL
  }else{
    stopifnot(length(dim(nodeData))==2)
    if(length(rownames(nodeData))==0){
      stop('nodeData must have rownames ',
      'and the rownames should be the names ofthe nodes.')
    }
    if(!all(node %in% rownames(nodeData))){
      warning('not all nodes is in the rownames of nodeData')
    }
  }
  preDefinedColor <- FALSE
  if(length(nodecolor) < 2){
    if(nodecolor %in% colnames(nodeData)){
      preDefinedColor <- TRUE
    }else{
      stop("nodecolor should have more than 1 elements")
    }
  }
  if(all(c('gene', 'miRNA') %in% colnames(cifNetwork))){
    if(length(setdiff(c('gene', 'miRNA'), names(nodeBorderColor))) > 0){
      stop("nodeBorderColor's element must be 'gene' and 'miRNA'")
    }
  }else{
    if(!'gene' %in% names(nodeBorderColor)){
      nodeBorderColor[["gene"]] <- 'darkgreen'
    }
  }
  
  if(length(edgeWeight)==1 && edgeWeight %in% colnames(cifNetwork)){
    cifNetwork[is.na(cifNetwork[, edgeWeight]), edgeWeight] <- 0
  }
  edL<-split(cifNetwork, cifNetwork[,"from"])
  edL<-lapply(node,function(.ele){
    .ele<-edL[[.ele]]
    if(is.null(.ele)){
      .ele<-list(edges=c(),weights=c())
    }else{
      if(length(edgeWeight)==1 && edgeWeight %in% colnames(.ele)){
        .ele<-list(edges=as.character(.ele$to),
                   weights=abs(.ele[, edgeWeight, drop=TRUE]))
      }else{
        .ele<-list(edges=as.character(.ele$to),weights=rep(1,length(.ele$to)))
      }
    }
  })
  names(edL)<-node
  gR<-new("graphNEL", nodes=node, edgeL=edL, edgemode="directed")
  ## set node default data
  nodeDataDefaults(gR, attr="label") <- NA
  nodeDataDefaults(gR, attr="fill")<-nodeBg
  nodeDataDefaults(gR, attr="size")<-nodesDefaultSize
  nodeDataDefaults(gR, attr='borderColor') <- nodeBorderColor[["gene"]]
  for(i in node) {
    nodeData(gR, n=i, attr="label") <- i
  }
  ## add additional message 
  additionalInfoColAll <- colnames(nodeData)
  additionalInfoColAll <-
    additionalInfoColAll[!additionalInfoColAll %in%
                           c("to", "from", "dir")]
  addAdditionalInfo <- function(gR, types, defaultValue){
    additionalInfoCol <- additionalInfoColAll[
      vapply(additionalInfoColAll, FUN=function(.e){
        inherits(nodeData[, .e], types) &&
          length(unique(nodeData[, .e])) > 1
      }, FUN.VALUE=FALSE)
    ]
    if(length(additionalInfoCol)){
      for(j in additionalInfoCol){
        nodeDataDefaults(gR, attr=j)<-defaultValue
        for(i in intersect(node, rownames(nodeData))){
          nodeData(gR, n=i, attr=j) <- nodeData[i, j]
        }
      }
    }
    return(gR)
  }
  ## additional characters
  gR <- addAdditionalInfo(gR, c("character", "factor"), "")
  ## additional numeric logical
  gR <- addAdditionalInfo(gR, c("numeric", "logical"), NA)
  ## set node size
  if(!'size' %in% colnames(nodeData)){
    for(i in unique(as.character(cifNetwork$from))){
      if(!is.na(i)){
        nodeData(gR, n=i, attr="size")<-
          ceiling(5*length(edL[[i]]$edges)/length(node)) *
          nodesDefaultSize/2 + nodesDefaultSize
      }
    }
  }
  
  ## set node color
  if(!'fill' %in% colnames(nodeData)){
    if(!preDefinedColor){
      if('logFC' %in% colnames(nodeData)){
        lfc <- nodeData[, 'logFC', drop=TRUE]
        lfcMax<-ceiling(max(abs(lfc[!is.infinite(lfc)]), na.rm = TRUE))
        lfcSeq<-seq(-1*lfcMax,lfcMax,length.out=length(nodecolor)+1)
        lfc[is.infinite(lfc)] <- sign(lfc) * lfcMax
        colors <- nodecolor[findInterval(lfc, lfcSeq, all.inside = TRUE)]
        names(colors) <- rownames(nodeData)
        colors[is.na(colors)] <- nodeBg
        for(i in intersect(node, names(colors))){
          nodeData(gR, n=i, attr="fill") <- colors[i]
        }
      }
    }else{
      colors <- nodeData[, nodecolor, drop=FALSE]
      for(i in intersect(node, rownames(colors))) {
        nodeData(gR, n=i, attr="fill") <- colors[i, nodecolor, drop=TRUE]
      }
    }
  }
  
  ## set node border color
  if('miRNA' %in% colnames(nodeData) &&
     'miRNA' %in% names(nodeBorderColor) &&
     !'borderColor' %in% colnames(nodeData)){
    if(is.logical(nodeData$miRNA)){
      miRNAs<-rownames(nodeData[nodeData[,"miRNA"], , drop=FALSE])
      for(i in intersect(node, miRNAs)){
        nodeData(gR, n=i, attr="borderColor")<-nodeBorderColor$miRNA
      }
    }
  }
  graph::nodeRenderInfo(gR) <- list(col=nodeData(gR, attr='borderColor'),
                                    fill=nodeData(gR, attr='fill'),
                                    ...)
  cifNetwork.s <- cifNetwork[!is.na(cifNetwork$from) & !is.na(cifNetwork$to),
                             , drop=FALSE]
  if(length(edgeWeight)==1 && edgeWeight %in% colnames(cifNetwork)){
    graph::edgeRenderInfo(gR) <- list(lwd=edgelwd)
    lwdScore <- cifNetwork.s[, edgeWeight, drop=TRUE]
    rg <- range(lwdScore)
    lwd <- findInterval(lwdScore, seq(from=rg[1], to=rg[2], length.out=10),
                        all.inside = TRUE)
    lwd <- lwd*edgelwd
    names(lwd) <- paste(cifNetwork.s$from, cifNetwork.s$to, sep='~')
    lwd <- lwd[names(graph::edgeRenderInfo(gR, 'lwd'))]
    graph::edgeRenderInfo(gR) <- list(lwd=lwd)
  }else{
    graph::edgeRenderInfo(gR) <- list(lwd=edgelwd)
  }
  
  gR
}