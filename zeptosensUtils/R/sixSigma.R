#' Peak selection w/ 6 sigma of base
#' 
#' @param arrayData FIXME
#' @param baseTreatment FIXME
#' @param nsamp number of samples 
#' @param antibodyNum number of antibodies used
#' @param sixsigname FIXME
#' 
#' @return FIXME
#' 
#' @concept zeptosens
#' @export
sixSigma <-
  function(arrayData,
           baseTreatment,
           nsamp,
           antibodyNum,
           sixsigname) {
    arrayData$readout <- log2(abs(arrayData$readout))
    
    base_read <- matrix(0, nrow = length(antibodyNum), ncol = nsamp)
    nnorm <- matrix(0, nrow = length(antibodyNum), ncol = 1)
    base_6sigma <- matrix(0, nrow = length(antibodyNum), ncol = 1)
    maxrdt <-  matrix(0, nrow = length(antibodyNum), ncol = 1)
    print(length(antibodyNum))
    col6sigma <- paste("antibody", "max_readout", "6sigma_dmso", "treatment", sep = "\t")
    
    write.table(
      col6sigma,
      file = sixsigname,
      row.names = F,
      col.names = F,
      sep = "\t",
      quote = F
    )
    
    for (i in 1:length(antibodyNum)) {
      for (j in 1:nsamp) {
        if (arrayData[(i - 1) * nsamp + j, "treatment"] == baseTreatment) {
          nnorm[i, 1] <- nnorm[i, 1] + 1
          
          base_read[i, nnorm[i, 1]] <- arrayData[(i - 1) * nsamp + j, "readout"]
          
        }
        #           print(arrayData[(i-1)*nsamp+j,])
      }
      base_6sigma[i, 1] <- 6 * sd(base_read[i, ], na.rm = T)
      t1 <- (i - 1) * nsamp + 1
      t2 <- i * nsamp
      array1 <- arrayData[t1:t2, ]
      maxrdt[i, 1] <- max(abs(array1[, "readout"]))
      
      print(paste(maxrdt[i, 1], base_6sigma[i, 1]))
      #        print((i-1)*nsamp+1)
      #        print(i*nsamp)
      if (maxrdt[i, 1] > base_6sigma[i, 1]) {
        write.table(
          paste(arrayData[(i - 1) * nsamp + 1, "antibody"], maxrdt[i, 1], base_6sigma[i, 1], sep = "\t"),
          file = sixsigname,
          append = T,
          quote = F,
          row.names = F,
          col.names = F
        )
      }
    }
  }
