#' Reads the Zeptosens File (.xls) into a Data.Frame
#' 
#' @param inputFile 
#' 
#' @return a data.frame where columns have been classified as "character" or "numeric"
#' 
#' @concept zeptosensPkg
#' @export
readZeptosensXls <- function(inputFile) {
    colClasses <- c("character", "numeric", "numeric", "character", "character", "character", 
                    "numeric", "character", "character", "character", "character", "character", 
                    "numeric", "numeric", "numeric", "character", "numeric", "numeric", 
                    "numeric", "numeric", "numeric", "numeric", "numeric", "character", 
                    "character", "character", "numeric")

    data <- read.xlsx2(inputFile, sheetIndex=1, colClasses=colClasses, stringsAsFactors=FALSE)
    
    return(data)
}