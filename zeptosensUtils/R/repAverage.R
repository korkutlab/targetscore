#' average replicates Zeptosens Data
#' 
#' @param array1 FIXME dataframe from readZeptosensExport
#' @param average_rep 
#' @param controlProbeIndex row number that that contains the first control probe
#' @param antibodyNum number of antibodies on array
#' 
#' @concept zeptosens
#' @export
repAverage <- function(array1) {
    
    array_rep1 <- array1[which(array1$replicate == "rep1"), ]
    array_rep2 <- array1[which(array1$replicate == "rep2"), ]
    array_rep3 <- array1[which(array1$replicate == "rep3"), ]
    array_rep4 <- array1[which(array1$replicate == "rep4"), ]
    average_rep <- array_rep1
    average_rep$readout <- (array_rep1$readout + array_rep2$readout + array_rep3$readout + 
        array_rep4$readout)/4
    nall <- nrow(array1)
    nuniq <- nall * 0.25
    
    for (i in 1:nuniq) {
        average_rep[i, ]$cv <- sd(c(array_rep1[i, ]$readout, array_rep2[i, ]$readout, 
            array_rep3[i, ]$readout, array_rep4[i, ]$readout))
        
    }
    
    # stdev <-
    # sd(c(array_rep1$readout,array_rep2$readout,array_rep3$readout,array_rep4$readout))
    # average_rep$cv <- stdev
    average_rep$replicate <- "average"
    
    return(average_rep)
} 
