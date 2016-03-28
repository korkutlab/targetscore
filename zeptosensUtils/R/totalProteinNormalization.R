#' Normalize Zeptosens Data wrt ttl protein level
#' 
#' @param array1 dataframe from readZeptosensExport function 
#' @param totalProteinLevels1 total protein levels for array 1 (use readZeptosensData format)
#' @param nsamples samples/Ab
#' @param antibodyNum number of antibodies on array
#' 
#' @details Normalize wrt ttl prt lvl
#' 
#' @concept zeptosens
#' @export
totalProteinNormalization <- function(Array1, totalProteinLevels1, nsamples, antibodyNum) {
    # tmpDf <- NULL
    normArray1 <- Array1
    for (i in 1:length(antibodyNum)) {
        for (j in 1:nsamples) {
            normArray1[(i - 1) * nsamples + j, "readout"] <- Array1[(i - 1) * nsamples + 
                j, "readout"]/totalProteinLevels1[j, 3]
        }
    }
    
    # normArray1 <- cbind(array1, normReadout=tmpDf[,'array1']) normArray2 <-
    # cbind(array2, normReadout=tmpDf[,'array2'])
    
    # tmp <- list(normArray1=normArray1) return(tmp)
    return(normArray1)
} 
