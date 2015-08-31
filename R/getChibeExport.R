genChibeExport <- function(antibodyMap=NULL) {
    if(is.null(antibodyMap)) {
        antibodyMap <- system.file("extdata", "antibodyMap.txt", package="zeptosensPkg")
    }
    
    read.table(antibodyMap, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
}