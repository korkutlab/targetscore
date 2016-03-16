#' Convert Intermediate Format to Matrix 
#' 
#' @param arrayData a data.frame as described in readZeptosensExport
#' @param colsForRowNames a vector of column numbers to be used for creating row names 
#' 
#' @return a matrix with antibodies as the column names and row names using the user selected columns 
#'   separated by underscores 
#'   
#' @concept zeptosensPkg  
#' @export
convertIntermediateFormatToMatrix <- function(arrayData, colsForRowNames) {
    mat <- matrix(NA, nrow=length(unique(arrayData$sNumber)), ncol=length(unique(arrayData$antibody)))
    for(i in 1:length(unique(arrayData$sNumber))) {
        cursNumber <- unique(arrayData$sNumber)[i]
        mat[i, ] <- arrayData[which(arrayData$sNumber == cursNumber), "readout"]
    }
    
    rowNames <- apply(arrayData[1:length(unique(arrayData$sNumber)), colsForRowNames], 1, function(x) {
        paste(x, collapse="_")    
    })
    
    colnames(mat) <- unique(arrayData$antibody)
    rownames(mat) <- rowNames
    
    return(mat) 
}
