#' Parse Causality Results to igraph
#'
#' @param results generateCausalityGraph() results; list with sif and format data.frames
#'
#' @return an igraph object
#'
#' @concept zeptosensPkg
#' @export
parseCausalityToIgraph <- function(results) {
  sif <- results$sif

  t1 <- sif[, c("PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B")]
  t2 <- t1[t1$INTERACTION_TYPE != "", ]
  g <- loadSifInIgraph(t2)

  # Shorten labels that have GBL
  labels <- V(g)$name
  idx <- which(grepl("_GBL", labels))
  labels[idx] <- strapplyc(V(g)$name[idx], "^(.*)_GBL.*$", simplify = TRUE)
  g <- set_vertex_attr(g, "label", index = V(g), value = labels)
  # g <- set_vertex_attr(g, "label", index = V(g), value=V(g)$name)

  # Get node colors
  formatTmp <- subset(results$format, componentProperty == "color" & componentType == "node", select = c(componentLabel, rgbColor))

  nodeColors <- NULL

  for (i in 1:length(V(g))) {
    idx <- which(formatTmp[, "componentLabel"] == V(g)$name[i])

    if (length(idx) == 1) {
      rgbTmp <- strsplit(formatTmp[idx, "rgbColor"], " ")[[1]]

      # Values must be between 0 and 1
      t1 <- as.numeric(rgbTmp) / 255
      hexColor <- do.call("rgb", as.list(t1))

      nodeColors <- c(nodeColors, hexColor)
    } else {
      nodeColors <- c(nodeColors, "#FFFFFF")
    }
  }

  # Get edge colors
  # formatTmp <- subset(results$format, componentProperty == "color" & componentType == "edge", select = c(componentLabel, rgbColor))

  edgeColors <- NULL
  edgeLty <- NULL

  for (i in 1:length(E(g))) {
    label <- paste(ends(g, E(g)[i])[1], E(g)[i]$interactionType, ends(g, E(g)[i])[2], sep = " ")
    # idx <- which(formatTmp[, "componentLabel"] == label)

    # if(length(idx) == 1) {
    # rgbTmp <- strsplit(formatTmp[idx, "rgbColor"], " ")[[1]]
    # t1 <- as.numeric(rgbTmp) / 255
    # hexColor <- do.call("rgb", as.list(t1))

    if (grepl("expression", E(g)[i]$interactionType)) {
      edgeLty <- c(edgeLty, 2)
    } else {
      edgeLty <- c(edgeLty, 1)
    }

    if (grepl("^phosphorylates", E(g)[i]$interactionType) || grepl("^upregulates", E(g)[i]$interactionType)) {
      edgeColors <- c(edgeColors, "#00ff00")
    } else {
      edgeColors <- c(edgeColors, "#ff0000")
    }
    # } else {
    #    edgeColors <- c(edgeColors, "#FFFFFF")
    # }
  }

  V(g)$color <- nodeColors
  E(g)$color <- edgeColors
  E(g)$lty <- edgeLty

  return(g)
}
