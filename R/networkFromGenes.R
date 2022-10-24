#' Build network by a list of given genes
#' @description By providing a list of given genes, build a network for input of filterNetwork.
#' @param genes A vector of character for interested genes.
#' @param interactionmap Transcription regulatory map. 
#'                       Column names of interactionmap must be 'from','to'
#' @param level Depth of node path
#' @param unrooted Return unrooted regulatory network table or not.
#' @return a list with elements:
#'  rootgene: The nodes with maximal connections.
#'  sifNetwork: Transcription regulatory network table.
#' @keywords network
#' @examples 
#' data("ce.interactionmap")
#' data("example.data")
#' genes <- as.character(example.data$ce.bind$from)
#' xx<-networkFromGenes(example.data$ce.bind, ce.interactionmap, level=2)
#' @export

networkFromGenes <- function(genes, interactionmap, level=3, unrooted=FALSE){
  stopifnot(!missing(genes))
  checkMCName(interactionmap)
  if(level>0){
    y<-interactionmap[interactionmap[ , "from"] %in%
                        unique(as.character(genes)),
                      1:2,drop=FALSE]
    TFbindingTable<-unique(y)
    level<-level-1
    if(level>0){
      y<-buildNetwork(TFbindingTable, interactionmap, level)
      TFbindingTable<-rbind(TFbindingTable, y)
      TFbindingTable<-unique(TFbindingTable)
    }
  }
  if(unrooted){
    return(list(rootgene=NULL, sifNetwork=TFbindingTable))
  }
  
  getParentNodes <- function(tfb, minParentNum){
    maxNode <- table(c(as.character(tfb$from), as.character(tfb$to)))
    maxNode <- sort(maxNode, decreasing = TRUE)
    maxNodeParents <- tfb[tfb$to %in% names(maxNode)[maxNode==max(maxNode)], ]
    parentNum <- length(unique(maxNodeParents$from))
    if(nrow(maxNodeParents)>1 && parentNum<minParentNum){
      maxNodeParents <- getParentNodes(maxNodeParents, parentNum)
    }
    maxNodeParents
  }
  getParentWithMaxChildren <- function(froms, tfb){
    tfb <- tfb[tfb$from %in% froms, , drop=FALSE]
    n <- split(tfb$to, tfb$from)
    n <- lapply(n, unique)
    names(n)[which.max(lengths(n))]
  }
  
  if(nrow(TFbindingTable)>0){
    p <- getParentNodes(TFbindingTable, 
                        minParentNum=length(unique(TFbindingTable$from)))
    p <- getParentWithMaxChildren(p$from, TFbindingTable)[1]
    list(rootgene=p, sifNetwork=TFbindingTable)
  }else{
    message("No transcription regulatory network.")
    list(rootgene=NULL, sifNetwork=NULL)
  }
}
