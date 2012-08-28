checkCName<-function(colName,exprsData){
    re<-TRUE
    if(is.numeric(colName)){
        if(colName>ncol(exprsData)){
            re<-FALSE
        }
    }else{
        if(!(colName %in% colnames(exprsData))){
            re<-FALSE
        }
    }
    re
}
checkMCName<-function(arr){
    cname<-c("from","to")
    if(!all(colnames(arr)[1:2]==cname)){
        stop(paste("colnames of",  deparse(substitute(arr)) , "must be c('from','to',...)"))
    }
}

checkMap<-function(interactionmap,x){
    cname<-c("from","to")
    checkMCName(interactionmap)
    checkMCName(x)
    tf_names<-unique(c(as.vector(t(interactionmap[,cname])),as.vector(x)))
    if(!all(!is.na(tf_names))){
        stop("NA is involved in input data.")
    }
}

checkmiRNAmap<-function(miRNAmap){
    if(dim(miRNAmap)[2] != 2) stop("miRNA map must be a two column matrix or data frame")
}

getMax<-function(dl, colname="logFC"){
    dl[match(max(abs(dl[ , colname])), abs(dl[ , colname])),]
}

getMedian<-function(dl, colname="logFC"){
    m<-median(dl[ , colname])
    dl[match(min(abs(m - dl[ , colname])), abs(m - dl[ , colname])), ]
}

getMin<-function(dl, colname="logFC"){
    dl[match(min(abs(dl[ , colname])), abs(dl[ , colname])), ]
}

inList<-function(needle,arr){
    toupper(needle) %in% toupper(arr)
}