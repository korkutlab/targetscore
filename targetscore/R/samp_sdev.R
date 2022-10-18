#' Calculate the standard deviation over all samples
#'
#' @param n_x vertically merged data from all samples rows condition columns proteins
#' @param n_sample samples number of input data.
#' @param n_prot Antibody number of input data.
#' @param n_dose doses number of input data.
#' @param replace_missing Replace missing values NA with 0 (DEFAULT: TRUE)
#' 
#' @examples 
#' # read proteomic response file
#' file <- system.file("test_data", "BT474.csv", package = "targetscore")
#' proteomic_responses <- read.csv(file, row.names = 1)
#' samp_d <- samp_sdev(
#' n_x = proteomic_responses, n_sample = dim(proteomic_responses)[1],
#' n_prot = dim(proteomic_responses)[2], n_dose = 1
#' )
#' 
#' @concept targetscore
#' @export
samp_sdev <- function(n_sample, n_prot, n_dose, n_x, replace_missing = TRUE) {
  n_x[is.na(n_x)] <- 0

  #    n_prot
  #    n_sample
  #    n_condition
  s_sdev <- array(0:0, dim = c(n_dose, n_prot))
  s_sdev2 <- array(0:0, dim = c(n_prot))

  #    for (i in 1:n_sample) {

  for (j in 1:n_prot) {
    for (k in 1:n_dose) {
      s_sdev[k, j] <- sd(n_x[((k - 1) * n_sample + 1):(k * n_sample), j])
    }
  }
  #    }
  #    for (i in 1:n_sample*n_dose) {

  for (j in 1:n_prot) {
    s_sdev2[j] <- sd(n_x[1:(n_sample * n_dose), j])
  }
  #    }
  return(s_sdev2)
}
