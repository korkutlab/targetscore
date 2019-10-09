#' average replicates Zeptosens Data
#'
#' @param array1 FIXME dataframe from readZeptosensExport
#' @param averageRep FIXME
#' @param controlProbeIndex row number that that contains the first control probe
#' @param antibodyNum number of antibodies on array
#'
#' @concept zeptosens
#' @export
repAverage <- function(array1) {
  arrayRep1 <- array1[which(array1$replicate == "rep1"), ]
  arrayRep2 <- array1[which(array1$replicate == "rep2"), ]
  arrayRep3 <- array1[which(array1$replicate == "rep3"), ]
  arrayRep4 <- array1[which(array1$replicate == "rep4"), ]
  averageRep <- arrayRep1
  averageRep$readout <- (arrayRep1$readout + arrayRep2$readout + arrayRep3$readout +
    arrayRep4$readout) / 4
  nall <- nrow(array1)
  nuniq <- nall * 0.25

  for (i in 1:nuniq) {
    averageRep[i, ]$cv <- sd(c(
      arrayRep1[i, ]$readout, arrayRep2[i, ]$readout,
      arrayRep3[i, ]$readout, arrayRep4[i, ]$readout
    ))
  }

  # stdev <-
  # sd(c(arrayRep1$readout,arrayRep2$readout,arrayRep3$readout,arrayRep4$readout))
  # averageRep$cv <- stdev
  averageRep$replicate <- "average"

  return(averageRep)
}
