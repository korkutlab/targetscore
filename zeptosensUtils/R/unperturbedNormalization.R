#' Normalize Zeptosens Data wrt unperturbed (or any other condition)
#' 
#' @param array1 data frame from readZeptosensExport function
#' @param normalFactor (treatment_name)
#' @param nSamples samples/Ab
#' @param antibodyNum number of antibodies on array
#' 
#' @details Normalize wrt ttl prt lvl
#' 
#' @concept zeptosens
#' @export
unperturbedNormalization <- function(array1, normalFactor, nSamples, antibodyNum) {
    print(normalFactor)
    # tmpDf <- NULL
    
    normArray1 <- array1
    norm_coef <- matrix(0, ncol = 1, nrow = antibodyNum)
    nnorm <- matrix(0, ncol = 1, nrow = antibodyNum)
    
    array1[is.na(array1[, "treatment"]), "treatment"] <- "NA"
    for (i in 1:antibodyNum) {
        for (j in 1:nSamples) {
            
            if (array1[(i - 1) * nSamples + j, "treatment"] == normalFactor) {
                
                norm_coef[i, 1] <- norm_coef[i, 1] + array1[(i - 1) * nSamples + 
                  j, "readout"]
                nnorm[i, 1] <- nnorm[i, 1] + 1
                print(nnorm[i, 1])
            }
        }
        
        norm_coef[i, 1] <- norm_coef[i, 1]/nnorm[i, 1]
    }
    for (i in 1:antibodyNum) {
        for (j in 1:nSamples) {
            normArray1[(i - 1) * nSamples + j, "readout"] <- array1[(i - 1) * nSamples + 
                j, "readout"]/norm_coef[i, 1]
        }
    }
    
    # normArray1 <- cbind(array1, normReadout=tmpDf[,'array1']) normArray2 <-
    # cbind(array2, normReadout=tmpDf[,'array2'])
    
    # tmp <- list(normArray1=normArray1) return(tmp)
    return(normArray1)
} 
