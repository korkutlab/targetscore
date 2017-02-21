#' Normalize Zeptosens Data
#' 
#' @param array1 FIXME
#' @param array2 FIXME
#' @param totalProteinLevels1 total protein levels for array 1 (use readZeptosensData format)
#' @param totalProteinLevels2 total protein levels for array 2 (use readZeptosensData format)
#' @param controlProbeIndex row number that that contains the first control probe
#' @param antibodyNum number of antibodies on array
#' 
#' @details Normalize data using inter-array method
#' 
#' @concept zeptosensPkg
#' @export
normalizeZeptosensData <- function(array1, array2, totalProteinLevels1, totalProteinLevels2, controlProbeIndex, antibodyNum) {
    ratios <- NULL
    
    for(i in 1:length(antibodyNum)) {
        curControlProbeIndex <- (controlProbeIndex*i)-1
        
        t1 <- sum(array1[c(curControlProbeIndex, curControlProbeIndex+1), "readout"])
        t2 <- sum(array2[c(curControlProbeIndex, curControlProbeIndex+1), "readout"])
        
        ratios <- c(ratios, t1/t2)
    }
    
    ratios <- data.frame(antibodyNames=antibodyNum, ratios=ratios)
    
    #write.table(ratios, file="inst/dataInst/a2058_skmel133_ctrlRatios.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    normCoef1 <- (ratios[, "ratios"]/ratios[, "ratios"]) %*% t(totalProteinLevels1[, 3])
    normCoef2 <- (ratios[, "ratios"]/ratios[, "ratios"]) %*% t(totalProteinLevels2[, 3])
    
    tmpDf <- NULL
    
    for(i in 1:nrow(normCoef1)) {
        for(j in 1:ncol(normCoef1)) {
            t1 <- array1[(i-1)*ncol(normCoef1)+j, "readout"] / normCoef1[i,j]
            t2 <- array2[(i-1)*ncol(normCoef1)+j, "readout"] / normCoef2[i,j]
            
            tmpDf <- rbind(tmpDf, data.frame(array1=t1, array2=t2))
        }
    }
    
    normArray1 <- cbind(array1, normReadout=tmpDf[,"array1"]) 
    normArray2 <- cbind(array2, normReadout=tmpDf[,"array2"]) 
    
    tmp <- list(normArray1=normArray1, normArray2=normArray2)
    return(tmp)
}