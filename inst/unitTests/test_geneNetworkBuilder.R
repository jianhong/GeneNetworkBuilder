test_uniqueExprsData<-function(){
    exprsData<-data.frame(symbols=c("A","B","C","A","B","A"),
                          FC=c(0,1,2,1,2,2))
    exprsDataMax<-uniqueExprsData(exprsData, "Max", "FC")
    exprsDataMed<-uniqueExprsData(exprsData, "Median", "FC")
    exprsDataMin<-uniqueExprsData(exprsData, "Min", "FC")
    checkEquals(nrow(exprsDataMax),3)
    checkEquals(exprsDataMax[exprsDataMax[,"symbols"]=="A", "FC"], 2)
    checkEquals(nrow(exprsDataMed),3)
    checkEquals(exprsDataMed[exprsDataMed[,"symbols"]=="A", "FC"], 1)
    checkEquals(nrow(exprsDataMin),3)
    checkEquals(exprsDataMin[exprsDataMin[,"symbols"]=="A", "FC"], 0)
}

test_convertID<-function(){
    idm<-c("A","B","C")
    names(idm)<-c("apple","boy","cap")
    x<-matrix(c("apple","Apple"),ncol=1,nrow=2)
    colnames(x)<-"ID"
    x<-convertID(x,idm,"ID")
    checkEquals(as.character(x[1,"ID"]),"A")
    checkTrue(is.na(x[2,"ID"]))
}

test_buildNetwork<-function(){
    interactionmap<-matrix(c("A","B",
                             "A","C",
                             "B","C",
                             "B","D",
                             "D","E",
                             "C","F",
                             "C","E"),
                           ncol=2,byrow=TRUE)
    colnames(interactionmap)<-c("from","to")
    x<-matrix(c("A","A"),ncol=2)
    colnames(x)<-c("from","to")
    xx<-buildNetwork(x,interactionmap,level=2)
    xxx<-buildNetwork(x,interactionmap,level=3)
    checkEquals(nrow(xx),7)
    checkEquals(nrow(xxx),8)
}

test_filterNetwork<-function(){
    rootgene<-"A"
    interactionmap<-matrix(c("A","B",
                             "A","C",
                             "B","C",
                             "B","D",
                             "C","E",
                             "C","F",
							 "C","H",
                             "D","E",
							 "D","I",
							 "F","G",
							 "F","K",
							 "G","L",
							 "I","J"),
                           ncol=2,byrow=TRUE)
    colnames(interactionmap)<-c("from","to")
    x<-matrix(c("A","A"),ncol=2)
    colnames(x)<-c("from","to")
    xx<-buildNetwork(x,interactionmap,level=5)
    exprsData<-data.frame(symbols=c(LETTERS[1:12]),
                          logFC=c(1,1,1,1,1,-1,0.2,-1,0.2,-1,1,1),
						  P.Value=c(0.001,0.001,0.001,0.001,0.001,0.1,0.001,0.001,0.001,0.001,0.001,0.001))
    miRNAlist<-c("C","I")
    if(.Platform$OS.type != "windows"){
        xxxf<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=FALSE,tolerance=0)
        xxxt<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=TRUE,tolerance=0)
        xxf<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=FALSE,tolerance=1)
        xxt<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=TRUE,tolerance=1)
        xf<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=FALSE,tolerance=2)
        xt<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=TRUE,tolerance=2)
        checkEquals(length(unique(unlist(xxxf[,1:2]))),8)
        checkEquals(length(unique(unlist(xxxt[,1:2]))),6)
        checkEquals(length(unique(unlist(xxf[,1:2]))),10)
        checkEquals(length(unique(unlist(xxt[,1:2]))),10)
        checkEquals(length(unique(unlist(xf[,1:2]))),12)
        checkEquals(length(unique(unlist(xt[,1:2]))),12)
        ## check miRNAlist could be missing
        null <- filterNetwork(rootgene,xx,exprsData,"symbols",tolerance=1)
    }
}

test_polishNetwork<-function(){
    rootgene<-"A"
    interactionmap<-matrix(c("A","B",
                             "A","C",
                             "B","C",
                             "B","D",
                             "D","E",
                             "C","E",
                             "C","F",
							 "C","H",
							 "F","G",
							 "D","I",
							 "I","J"),
                           ncol=2,byrow=TRUE)
    colnames(interactionmap)<-c("from","to")
    x<-matrix(c("A","A"),ncol=2)
    colnames(x)<-c("from","to")
    xx<-buildNetwork(x,interactionmap,level=3)
    exprsData<-data.frame(symbols=c(LETTERS[1:10]),
                          logFC=c(1,1,1,1,1,-1,0.2,-1,0.2,-1),
						  P.Value=c(0.001,0.001,0.001,0.001,0.001,0.1,0.001,0.001,0.001,0.001))
    miRNAlist<-c("C","I")
    
    if(.Platform$OS.type != "windows"){
        xxxf<-filterNetwork(rootgene,xx,exprsData,"symbols",miRNAlist,miRNAtol=FALSE)
        
        xxxf$info1 <- sample(c("groupA", "groupB"),
                             size = nrow(xxxf),
                             replace = TRUE)
        xxxf$info2 <- sample(c(FALSE, TRUE),
                             size = nrow(xxxf),
                             replace = TRUE)
        xxxf$info3 <- sample(seq.int(7),
                             size = nrow(xxxf),
                             replace = TRUE)
        xxxf$info4 <- factor(sample(LETTERS,
                                    size = nrow(xxxf),
                                    replace = TRUE))
        gR<-polishNetwork(xxxf)
        l<-RBGL::sp.between(gR,"A","E")
        checkEquals(l[[1]]$length,3)
    }
}