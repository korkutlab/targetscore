#' Match a vector of genes to a SIF network
#'
#' @param genes1 a vector of genes
#' @param genes2 an optional vector if the source
#' @param annot_edgelist a data.frame; the first two columns are interaction
#'   participants and the third column is an optional annotation
#' @param antibody_vec a vector with the names of the antibodies (e.g. colnames(proteomicResponses));
#'   indicies returned from this function will be mapped to this vector
#' @param use_annot a boolean as to whether the optional annotation column values
#'   in the annot_edgelist should outputted in the results
#' @param verbose a boolean to show debugging information
#'
#' @return a data.frame, columns 1-2 are indicies of the edgelist participants,
#'   column 3 is the annotation values for the edgelist, 4-5 the names of the
#'   edgelist participants
#'
#' @concept zeptosensPkg
#' @export
matchGenesToEdgelist <- function(genes1, genes2 = NULL, annot_edgelist, antibody_vec, 
                                 use_annot = FALSE, verbose = FALSE) {
  results <- data.frame(
    gene1 = numeric(), gene2 = numeric(), annot = numeric(),
    gene1Name = character(), gene2Name = character(),
    stringsAsFactors = FALSE
  )

  if (is.null(genes2)) {
    genes2 <- genes1
  }

  if (is.null(names(genes1)) || is.null(names(genes2))) {
    stop("ERROR: genes1 and genes2 must be named vectors with antibody names")
  }

  for (i in 1:length(genes1)) {
    # if(verbose) {
    #    cat("I: ", i, "\n")
    # }

    tmp_idx <- which(annot_edgelist[, 1] == genes1[i])
    cur_annot_edgelist <- annot_edgelist[tmp_idx, ]

    for (j in 1:length(genes2)) {
      # if(verbose) {
      #    cat("J: ", j, "\n")
      # }

      tmp_idx <- which(cur_annot_edgelist[, 1] == genes1[i] & cur_annot_edgelist[, 2] == genes2[j])

      if (length(tmp_idx) == 1) {
        if (use_annot) {
          annot <- cur_annot_edgelist[tmp_idx, 3]
        } else {
          annot <- NA
        }

        # Get indicies based off antibody names rather than the gene names
        gene1_ab_idx <- which(antibody_vec == names(genes1[i]))
        gene2_ab_idx <- which(antibody_vec == names(genes2[j]))

        tmp_results <- data.frame(
          gene1 = gene1_ab_idx,
          gene2 = gene2_ab_idx,
          annot = annot,
          gene1Name = genes1[i],
          gene2Name = genes2[j],
          gene1Ab = names(genes1[i]),
          gene2Ab = names(genes2[j]),
          stringsAsFactors = FALSE
        )
        results <- rbind(results, tmp_results)
      }

      if (length(tmp_idx) > 1) {
        stop("ERROR: Multiple shortest paths found. Check SignedPC network.")
      }
    }
  }

  return(results)
}
