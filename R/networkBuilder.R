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
        y<-interactionmap[interactionmap[ , "from"] %in% unique(as.character(TFbindingTable[ , "to"])), 1:2,drop=F]
        y<-unique(y)
        z<-y[!(y[,"to"] %in% TFbindingTable[,"to"]), , drop=F]
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
#' @param cifNetwork dataframe used to draw network graph. column names of 
#'                   cifNetwork must contain 'from', 'to', 'logFC' and 'miRNA'
#' @param nodesDefaultSize nodes default size
#' @param useLogFCAsWeight how to determine the weights for each nodes. 
#'                         If TURE, use logFC value as weight.
#'                         If FALSE, use constant 1 as weight.
#' @param nodecolor a character vector of color set. 
#'                  The node color will be mapped to color set by log fold change.
#'                  Or the column names for the colors.
#' @param nodeBg background of node
#' @param nodeBorderColor a list of broder node color set. 
#'                        nodeBorderColor's element must be gene and miRNA
#' @param edgelwd the width of edge
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
                        nodesDefaultSize=48, useLogFCAsWeight=FALSE, 
                        nodecolor=colorRampPalette(c("green", "yellow", "red"))(5), nodeBg="white",
                        nodeBorderColor=list(gene='darkgreen',miRNA='darkblue'), 
                        edgelwd=0.25, ...)
{
  cname<-c("from", "to")
  if(!is.data.frame(cifNetwork)){
    stop("cifNetwork should be a data.frame")
  }
  if(length(intersect(c("from", "to", "logFC", "miRNA"), colnames(cifNetwork)))<4){
    stop("colnames of cifNetwork must contain 'from', 'to', 'logFC' and 'miRNA'");
  }
  preDefinedColor <- FALSE
  if(length(nodecolor) < 2){
    if(nodecolor %in% colnames(cifNetwork)){
      preDefinedColor <- TRUE
    }else{
      stop("nodecolor should have more than 1 elements")
    }
  }
  if(length(setdiff(c('gene', 'miRNA'), names(nodeBorderColor))) > 0){
    stop("nodeBorderColor's element must be 'gene' and 'miRNA'")
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
  edL<-split(cifNetwork[,c("to","logFC")],cifNetwork[,"from"])
  edL<-lapply(node,function(.ele,edL,useLogFCAsWeight){
    .ele<-edL[[.ele]]
    if(is.null(.ele)){
      .ele<-list(edges=c(),weights=c())
    }else{
      if(useLogFCAsWeight){
        .ele<-list(edges=as.character(.ele$to),weights=abs(.ele$logFC))
      }else{
        .ele<-list(edges=as.character(.ele$to),weights=rep(1,length(.ele$to)))
      }
    }
  },edL,useLogFCAsWeight)
  names(edL)<-node
  gR<-new("graphNEL", nodes=node, edgeL=edL, edgemode="directed")
  ## set node default data
  nodeDataDefaults(gR, attr="label") <- NA
  nodeDataDefaults(gR, attr="logFC") <- 0
  nodeDataDefaults(gR, attr="miRNA") <- FALSE
  for(i in node) {
    nodeData(gR, n=i, attr="label") <- i
    nodeData(gR, n=i, attr="logFC") <-
      cifNetwork[match(i, cifNetwork$to), "logFC"]
    nodeData(gR, n=i, attr="miRNA") <-
      cifNetwork[match(i, cifNetwork$to), "miRNA"]
  }
  ## set node size
  nodeDataDefaults(gR, attr="size")<-nodesDefaultSize
  for(i in unique(as.character(cifNetwork$from))){
    if(!is.na(i)){
      nodeData(gR, n=i, attr="size")<-
        ceiling(5*length(edL[[i]]$edges)/length(node)) *
        nodesDefaultSize/2 + nodesDefaultSize
    }
  }
  ## add additional message 
  additionalInfoColAll <- colnames(cifNetwork)
  additionalInfoColAll <-
    additionalInfoColAll[!additionalInfoColAll %in%
                        c("to", "from", "gene", "P.Value", "logFC", "miRNA",
                          "dir")]
  addAdditionalInfo <- function(gR, types, defaultValue){
    additionalInfoCol <- additionalInfoColAll[
      vapply(additionalInfoColAll, FUN=function(.e){
        inherits(cifNetwork[, .e], types) &&
          length(unique(cifNetwork[, .e])) > 1
      }, FUN.VALUE=FALSE)
    ]
    if(length(additionalInfoCol)){
      for(j in additionalInfoCol){
        nodeDataDefaults(gR, attr=j)<-defaultValue
        for(i in unique(as.character(cifNetwork$from))){
          nodeData(gR, n=i, attr=j) <-
            cifNetwork[match(i, cifNetwork$to), j]
        }
      }
    }
    return(gR)
  }
  ## additional characters
  gR <- addAdditionalInfo(gR, c("character", "factor"), "")
  ## additional numeric logical
  gR <- addAdditionalInfo(gR, c("numeric", "logical"), NA)
  
  ## set node color    
  nodeDataDefaults(gR, attr="fill")<-nodeBg
  if(!preDefinedColor){
    lfcMax<-ceiling(max(abs(cifNetwork[!is.na(cifNetwork$logFC),"logFC"])))
    lfcSeq<-seq(-1*lfcMax,lfcMax,length.out=length(nodecolor)+1)
    colset<-unique(cifNetwork[!is.na(cifNetwork$logFC),c("to","logFC")])
    colset<-apply(colset, 1, function(.ele,color,lfcSeq){
      id=0
      for(i in 1:length(lfcSeq)){
        .elelfc<-as.numeric(as.character(.ele[2]))
        if(lfcSeq[i]<=.elelfc & lfcSeq[i+1]>=.elelfc){
          id=i
          break
        }
      }
      if(id!=0){
        c(.ele,nodecolor[id])
      }else{
        c(.ele,nodeBg)
      }
    },nodecolor,lfcSeq)
    colors<-colset[3,]
    names(colors)<-colset[1,]
    for(i in names(colors)){
      nodeData(gR, n=i, attr="fill")<-colors[i]
    }
    colset<-node[!node %in% names(colors)]
    names(colset)<-colset
    colset<-nodeBg
    colors<-c(colors,colset)
  }else{
    colors <- cifNetwork[match(node, cifNetwork$to), nodecolor]
    names(colors) <- node
    colors[is.na(colors)] <- nodeBg
    for(i in node) {
      tmp <- cifNetwork[match(i, cifNetwork$to), nodecolor]
      nodeData(gR, n=i, attr="fill") <- 
        ifelse(is.na(tmp), nodeBg, tmp)
    }
  }
  ## set node border color    
  miRNAs<-unique(as.character(cifNetwork[cifNetwork[,"miRNA"],"to"]))
  nodeBC<-character(length(node))
  names(nodeBC)<-node
  nodeDataDefaults(gR, attr="borderColor")<-nodeBorderColor$gene
  for(i in node) {
    if(i %in% miRNAs){
      nodeBC[i]<-nodeBorderColor$miRNA
      nodeData(gR, n=i, attr="borderColor")<-nodeBorderColor$miRNA
    }else{
      nodeBC[i]<-nodeBorderColor$gene
    }
  }
  
  graph::nodeRenderInfo(gR) <- list(col=nodeBC, fill=colors, ...)
  graph::edgeRenderInfo(gR) <- list(lwd=edgelwd)
  gR
}