#' Normalize Zeptosens Data wrt unperturbed (or any other condition)
#' 
#' @param array1 FIXME
#' @param  normal_factor (treatment_name)
#' @param nsamples samples/Ab
#' @param antibodyNum number of antibodies on array
#' 
#' @details Normalize wrt ttl prt lvl
#' 
#' @concept zeptosensPkg
#' @export
unperturbednormalization <- function(Array1, normal_factor, nsamples, antibodyNum) {
    print(normal_factor)
    tmpDf <- NULL
    
    normArray1 <- Array1
    norm_coef <- matrix(0,ncol=1,nrow=length(antibodyNum))
    nnorm <- matrix(0,ncol=1,nrow=length(antibodyNum))
    
    Array1[is.na(Array1[,"treatment"]), "treatment"] <- "NA"
    for(i in 1:length(antibodyNum)) {
        for(j in 1:nsamples) {
            
            if(Array1[(i-1)*nsamples+j, "treatment"]==normal_factor){

                norm_coef[i,1] = norm_coef[i,1]+Array1[(i-1)*nsamples+j, "readout"] 
                nnorm[i,1]=nnorm[i,1]+1
            }
        }
        norm_coef[i,1]=norm_coef[i,1]/nnorm[i,1]
    }    
    for(i in 1:length(antibodyNum)) {
        for(j in 1:nsamples) {
            normArray1[(i-1)*nsamples+j, "readout"] <- Array1[(i-1)*nsamples+j, "readout"] / norm_coef[i,1]
        }
    }
    
#    normArray1 <- cbind(array1, normReadout=tmpDf[,"array1"]) 
#    normArray2 <- cbind(array2, normReadout=tmpDf[,"array2"]) 
    
    tmp <- list(normArray1=normArray1)
#    return(tmp)
    return(normArray1)
}