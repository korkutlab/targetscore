#' Match a vector of genes to a SIF network
#' 
#' @param genes1 a vector of genes 
#' @param genes2 an optional vector if the source
#' @param annotEdgelist a data.frame; the first two columns are interaction 
#'   participants and the third column is an optional annotation
#' @param antibodyVec a vector with the names of the antibodies (e.g. colnames(proteomicResponses)); 
#'   indicies returned from this function will be mapped to this vector
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
matchGenesToEdgelist <- function(genes1, genes2=NULL, annotEdgelist, antibodyVec, useAnnot=FALSE, verbose=FALSE) {
    results <- data.frame(gene1=numeric(), gene2=numeric(), annot=numeric(), 
                              gene1Name=character(), gene2Name=character(), 
                              stringsAsFactors=FALSE)
    
    if(is.null(genes2)) {
        genes2 <- genes1
    } 
    
    if(is.null(names(genes1)) || is.null(names(genes2))) {
        stop("ERROR: genes1 and genes2 must be named vectors with antibody names")
    }
    
    for(i in 1:length(genes1)) {
        #if(verbose) {
        #    cat("I: ", i, "\n")
        #}
        
        tmpIdx <- which(annotEdgelist[,1] == genes1[i])
        curAnnotEdgelist <- annotEdgelist[tmpIdx, ]
        
        for(j in 1:length(genes2)) {
            #if(verbose) {
            #    cat("J: ", j, "\n")
            #}
            
            tmpIdx <- which(curAnnotEdgelist[,1] == genes1[i] & curAnnotEdgelist[,2] == genes2[j])
            
            if(length(tmpIdx) == 1) {
                if(useAnnot) {
                    annot <- curAnnotEdgelist[tmpIdx, 3]
                } else {
                    annot <- NA
                }
                
                # Get indicies based off antibody names rather than the gene names
                gene1AbIdx <- which(antibodyVec == names(genes1[i]))
                gene2AbIdx <- which(antibodyVec == names(genes2[j]))
                
                tmpResults <- data.frame(gene1=gene1AbIdx, 
                                         gene2=gene2AbIdx, 
                                         annot=annot, 
                                         gene1Name=genes1[i], 
                                         gene2Name=genes2[j], 
                                         gene1Ab=names(genes1[i]), 
                                         gene2Ab=names(genes2[j]),
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