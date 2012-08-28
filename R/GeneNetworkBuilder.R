uniqueExprsData<-function(exprsData, method, condenseName='logFC'){
    if(!(method %in% c("Max", "Median", "Min"))){
        stop("method must be Max, Median or Min")
    }
    if(!GeneNetworkBuilder:::checkCName("symbols", exprsData)){
        stop("symbols is not a valide colname of exprsData")
    }
    if(!GeneNetworkBuilder:::checkCName(condenseName, exprsData)){
        stop(paste(condenseName," is not a valide colname of exprsData"))
    }
    if(!is.data.frame(exprsData)){
        exprsData<-as.data.frame(exprsData)
    }
    if(!is.numeric(exprsData[ , condenseName])){
        stop(paste("class of", condenseName, "is not a numeric column"))
    }
    exprsData<-switch(method,
                           Max       =plyr::ddply(exprsData, plyr::.(symbols), GeneNetworkBuilder:::getMax, condenseName),
                           Median    =plyr::ddply(exprsData, plyr::.(symbols), GeneNetworkBuilder:::getMedian, condenseName),
                           Min       =plyr::ddply(exprsData, plyr::.(symbols), GeneNetworkBuilder:::getMin, condenseName)
                           )
    exprsData
}

convertID<-function(x, IDsMap, ByName=c("from", "to")){
    if((!is.character(IDsMap)) | (is.null(IDsMap))){
        stop("invalide IDsMap")
    }
    for(i in 1:length(ByName)){
        if(!GeneNetworkBuilder:::checkCName(ByName[i],x)){
            stop(paste(ByName[i],"is not a valide colname of x"))
        }
        x[,ByName[i]]<-IDsMap[as.character(x[,ByName[i]])]
    }
    x
}

buildNetwork<-function(TFbindingTable, interactionmap, level=3){
    GeneNetworkBuilder:::checkMap(interactionmap, TFbindingTable)
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

filterNetwork<-function(rootgene, sifNetwork, exprsData, mergeBy="symbols", miRNAlist, remove_miRNA=FALSE,
                    tolerance=0, cutoffPVal=0.01, cutoffLFC=0.5, minify=TRUE, miRNAtol=FALSE)
{
    GeneNetworkBuilder:::checkMCName(sifNetwork)
    if(!is.vector(miRNAlist)){
        stop("miRNAlist should be a vector")
    }
    if(!GeneNetworkBuilder:::checkCName(mergeBy, exprsData)){
        stop(paste(mergeBy, "is not a column name of exprsData"))
    }
    if(!GeneNetworkBuilder:::checkCName("logFC", exprsData)){
        stop("logFC is not a column name of exprsData")
    }
    if(!is.numeric(exprsData[ , "logFC"])){
        stop("class of exprsData[ , \"logFC\"] is not a numeric column")
    }
    if(!GeneNetworkBuilder:::checkCName("P.Value", exprsData)){
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
    cifNetwork<-merge(sifNetwork, exprsData, by.x="to", by.y=mergeBy, all.x=TRUE)
##   convert NA to 0 for logFC
    cifNetwork.logFC<-cifNetwork[,"logFC"]
    cifNetwork.logFC[is.na(cifNetwork.logFC)]<-0.0
    cifNetwork.pValue<-cifNetwork[,"P.Value"]
    cifNetwork.pValue[is.na(cifNetwork.pValue)]<-0.0
##   label microRNA
    cifNetwork$miRNA<-FALSE
    cifNetwork$dir<-2
    if(length(miRNAlist)>0){
        cifNetwork$miRNA<-ifelse(cifNetwork$to %in% miRNAlist, TRUE, FALSE)
        cifNetwork$dir<-ifelse(cifNetwork$from %in% miRNAlist, 0, 2)
    }
##   remove micorRNA
    if(remove_miRNA){
        cifNetwork<-cifNetwork[!cifNetwork$miRNA, ]
        cifNetwork.logFC<-cifNetwork.logFC[!cifNetwork$miRNA]
    }
    rootlogFC<-exprsData[exprsData[ , mergeBy] == rootgene, "logFC"]
    rootlogFC<-rootlogFC[!is.na(rootlogFC)]
    rootlogFC<-ifelse(length(rootlogFC) < 1, 0.0, rootlogFC[1])
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
    cifNetwork <- merge(cifNetwork, cifNetwork.list)
    cifNetwork
}

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
    if(length(nodecolor) < 2){
        stop("nodecolor should have more than 1 elements")
    }
    if(length(setdiff(c('gene', 'miRNA'), names(nodeBorderColor))) > 0){
        stop("nodeBorderColor's element must be 'gene' and 'miRNA'")
    }
    cifNetwork<-cifNetwork[!duplicated(cifNetwork[,cname]), ]
    edge<-cifNetwork[cifNetwork$from!="" & cifNetwork$to!="", cname]
    node<-c(as.character(unlist(edge)))
    node<-node[!is.na(node)]
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
## set node size
    nodeDataDefaults(gR, attr="size")<-nodesDefaultSize
    for(i in unique(as.character(cifNetwork$from))){
        nodeData(gR, n=i, attr="size")<-ceiling(5*length(edL[[i]]$edges)/length(node)) * nodesDefaultSize/2 + nodesDefaultSize
    }
## set node color    
    nodeDataDefaults(gR, attr="fill")<-nodeBg
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