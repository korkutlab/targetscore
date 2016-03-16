#' Generate ChiBE Export 
#'
#' @param antibodies a vector of strings of antibody names 
#' @param antibodyMap a string, antibody map file name (default: NULL, if NULL then use zentosensPkg internal antibody map)
#' @param chibeFilename a string, the name of the file with gene/site information needed for ChiBE
#' 
#' @concept zeptosensPkg
#' @export
genChibeExport <- function(antibodies, antibodyMapFilename=NULL, chibeFilename="abChibe.txt") {
    if(is.null(antibodyMapFilename)) {
        antibodyMapFilename <- system.file("extdata", "antibodyMap.txt", package="zeptosensPkg")
    }
    
    antibodyMap <- read.table(antibodyMapFilename, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    idx <- which(antibodies %in% antibodyMap[, "AntibodyLabel"])
    
    # Make lower case and split into synonyms 
    tmp <- tolower(antibodyMap[,"AntibodyLabel"])
    namesLst <- strsplit(tmp, "\\|")
    
    results <- list()
        
    for(curAntibody in antibodies) {
        q <- tolower(curAntibody)
        curResults <- searchListOfVectors(q, namesLst)
        results[[curAntibody]] <- as.vector(curResults)
    }
    
    results <- unlist(results)
    
    notFoundAntibodies <- names(which(is.na(results)))
    warning("Antibodies Not Found: ", paste(notFoundAntibodies, collapse=", "))

    tmp <- as.numeric(results)
    idx <- tmp[!is.na(tmp)]
    
    chibeData <- antibodyMap[idx, c("NodeName", "Symbol", "Sites", "Effect")]
    
    write.table(chibeData, file=chibeFilename, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    return(chibeData)
}

#' Search List of Vectors
#'
#' @param q query vector
#' @param lst list of vectors to search
#'
#' @return a list of vectors with the same length as the query vector, each list
#'   entry will have indicies for lst where there was a match with the query
#'   vector. Return NA if there were no matches.
#'
#' @details
#' Taken from: http://stackoverflow.com/questions/11002391/fast-way-of-getting-index-of-match-in-list
#'
#' @examples
#' lst <- list(1:3, 3:5, 3:7)
#' q <- c(3, 5)
#' results <- searchListOfVectors(q, lst)
#' names(results) <- q
#'
#' @concept rcellminerUtils
#' @export
searchListOfVectors <- function(q, lst) {
    tmp <- rep(seq_along(lst), sapply(lst, length))
    resultsSe <- sapply(q, function(x) tmp[which(unlist(lst) %in% x)])
    
    if(class(resultsSe) == "list") {
        return(NA)
    } else {
        return(resultsSe)
    }
}