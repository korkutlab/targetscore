#standard deviation over all samples at each condition
#' @param nX vertically merged data from all samples rows: condition columns: proteins
#' @param sSdev standard deviation of response over all samples for each protein in each condition 
#' @param nSample number of samples
#' @param nProt number of proteins
#' @param nCond number of consitions


sampSdev <- function(nX,sSdev,nSample,nProt,nCond)
    
#    nProt
#    nSample
#    nCondition        

    for (j in nProt){
        
        for (i in 1:nSample){
            for (k in 1:nCond){
                sSdev[i,j] <- sd(nX[((i-1)*nCond+1):i*nCond,j])
            }
    }
}
