#' Match a vector of genes to a SIF network
#' 
#' @param genes1 a vector of genes 
#' @param genes2 an optional vector if the source
#' @param annotEdgelist a data.frame; the first two columns are interaction 
#'   participants and the third column is an optional annotation
#' @param useAnnot a boolean as to whether the optional annotation column values 
#'   in the annotEdgelist should outputted in the results
#' @param verbose a boolean to show debugging information  
#'
#' @return a data.frame, columns 1-2 are indicies of the edgelist participants, 
#'   column 3 is the annotation values for the edgelist, 4-5 the names of the 
#'   edgelist participants
#'
#' @concept zeptosensPkg
#' @export
matchGenesToEdgelist <- function(genes1, genes2=NULL, annotEdgelist, useAnnot=FALSE, verbose=FALSE) {
    results <- data.frame(gene1=numeric(), gene2=numeric(), annot=numeric(), 
                              gene1Name=character(), gene2Name=character(), 
                              stringsAsFactors=FALSE)
    
    if(is.null(genes2)) {
        genes2 <- genes1
    } 
    
    for(i in 1:length(genes1)) {
        if(verbose) {
            cat("I: ", i, "\n")
        }
        
        for(j in 1:length(genes2)) {
            #cat("A")
            gene1 <- genes1[i]
            gene2 <- genes2[j]
            
            tmpIdx <- which(annotEdgelist[,1] == gene1 & annotEdgelist[,2] == gene2)
            
            if(length(tmpIdx) == 1) {
                if(useAnnot) {
                    annot <- annotEdgelist[tmpIdx, 3]
                } else {
                    annot <- NA
                }
                
                tmpResults <- data.frame(gene1=i, 
                                         gene2=j, 
                                         annot=annot, 
                                         gene1Name=gene1, 
                                         gene2Name=gene2, 
                                         stringsAsFactors=FALSE)
                results <- rbind(results, tmpResults)
            } 
            
            if (length(tmpIdx) > 1) {
                stop("ERROR: Multiple shortest paths found. Check SignedPC network.")
            }
        }
    }
    
    return(results)
}