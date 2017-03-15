#' Reads the Zeptosens File (.xls) into a Data.Frame
#' 
#' @param inputFile 
#' 
#' @return a data.frame where columns have been classified as 'character' or 'numeric'
#' 
#' @concept zeptosens
#' @export
readZeptosensXls <- function(inputFile) {
    colClasses <- c("character", "numeric", "numeric", "character", "character", 
        "character", "numeric", "character", "character", "character", "character", 
        "character", "numeric", "numeric", "numeric", "character", "numeric", "numeric", 
        "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character", 
        "character", "numeric")
    
    data <- tryCatch({
        data <- xlsx::read.xlsx2(inputFile, sheetIndex = 1, colClasses = colClasses, 
            stringsAsFactors = FALSE)
    }, error = function(err) {
        stop(paste("ERROR:  ", err, "\nTIP: The file may not be a valid XLS or XLSX file. Check that you cannot read the file in a text editor."))
    })
    
    rownames(data) <- 1:nrow(data)
    
    return(data)
} 
