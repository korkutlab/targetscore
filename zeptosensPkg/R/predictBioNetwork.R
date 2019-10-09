#' Predict Network from signedPC (Bio-inferred network)
#'
#' @param nProt TBA
#' @param proteomicResponses TBA
#' @param maxDist TBA (default: 1)
#' @param antibodyMapFile a listing of antibodies, their associated genes, and modification sites
#' @param distFile A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#'
#' @examples
#'
#' @concept zeptosensPkg
#' @export
predictBioNetwork <- function(nProt, proteomicResponses, maxDist, 
                              antibodyMapFile = NULL, distFile = NULL, verbose = F) {

  # match Ab names to gene names & posttranslational modifications
  if (is.null(antibodyMapFile)) {
    antibodyMapFile <- system.file("targetScoreData", "antibodyMap.txt", package = "zeptosensPkg")
  }
  mabToGenes <- antibodyMapFile

  if (verbose) {
    print(mabToGenes)
  }

  # pathway distance matrix
  if (is.null(distFile)) {
    distFile <- system.file("targetScoreData", "distances.txt", package = "zeptosensPkg")
  }

  tmpDist <- read.table(distFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  if (nProt != ncol(proteomicResponses)) {
    stop("ERROR: nProt is not equal to proteomicResponses column number")
  }

  # Filter dist to only keep those with a distance less than maxDist
  idx <- which(tmpDist[, 3] <= maxDist)
  dist <- tmpDist[idx, ]
  # dist
  idxAbMap <- which(mabToGenes[, 1] %in% colnames(proteomicResponses))
  #    print(mabToGenes[, 1])
  #    print(colnames(proteomicResponses))
  #    print(unique(mabToGenes[idxAbMap,]))
  #    if(length(idxAbMap) < nProt) {
  #    print(length(unique(mabToGenes[idxAbMap,1])))
  if (length(idxAbMap) < nProt) {
    #        print(length(idxAbMap))
    stop("ERROR: Not all columns in data were matched in antibody map")
  }
  #    print((unique(mabToGenes[idxAbMap, 1])))
  if (length(unique(mabToGenes[idxAbMap, 1])) != nProt) {
    print(unique(mabToGenes[idxAbMap, 1]))
    stop("ERROR: Mismatch in the number of selected antibodies and the number of proteomic responses")
  }

  # Used by matchGenesToEdgelist to account for the cases where multiple entries for the same antibody
  # exist in the antibody map
  antibodyMapSubset <- mabToGenes[idxAbMap, ]

  mabGenes <- mabToGenes[idxAbMap, 4]
  names(mabGenes) <- mabToGenes[idxAbMap, 1]

  distList <- matchGenesToEdgelist(
    genes1 = mabGenes, genes2 = NULL, annotEdgelist = dist,
    antibodyVec = colnames(proteomicResponses), useAnnot = TRUE, verbose = TRUE
  )


  # # Get interactions for measured genes
  # dist_gene1 <- pmatch(dist[, 1], mabToGenes[measured_genes, 4], duplicates.ok = TRUE)
  # dist_gene1Name <- mabToGenes[measured_genes, 4][dist_gene1]
  # dist_gene2 <- pmatch(dist[, 2], mabToGenes[measured_genes, 4], duplicates.ok = TRUE)
  # dist_gene2Name <- mabToGenes[measured_genes, 4][dist_gene2]
  #
  # # distance framework
  # distList <- data.frame(dist_gene1=dist_gene1, dist_gene2=dist_gene2, dist=dist[, 3],
  #                         dist_gene1Name=dist_gene1Name, dist_gene2Name=dist_gene2Name, stringsAsFactors = FALSE)
  # distList[is.na(distList[, 1]), 3] <- 100
  # distList[is.na(distList[, 2]), 3] <- 100
  # distList[, 1]
  
  # distInd(upstream,downstream)
  distInd <- matrix(Inf, ncol = nProt, nrow = nProt, 
                     dimnames = list(colnames(proteomicResponses), colnames(proteomicResponses))) 

  for (i in 1:length(distList[, 1])) {
    distInd[distList[i, 1], distList[i, 2]] <- distList[i, 3]

    if (distInd[distList[i, 1], distList[i, 2]] > maxDist) {
      distInd[distList[i, 1], distList[i, 2]] <- Inf
    }

    if (distInd[distList[i, 1], distList[i, 2]] == 0) {
      distInd[distList[i, 1], distList[i, 2]] <- Inf
    }
  }
  write.table(distInd, file = "distInd.txt", quote = F)
  # cov318 results ~2100

  ### get the network product### phospFile <- system.file('SignedPC', 'phosphorylates.txt',
  ### package='zeptosensPkg') phosp <- read.csv(phospFile, sep='\t', header=TRUE, na.strings =
  ### c('', ' ')) phosp <- phosp[,-3] dephospFile <- system.file('SignedPC',
  ### 'dephosphorylates.txt', package='zeptosensPkg') dephosp <- read.csv(dephospFile, sep='\t',
  ### header=TRUE, na.strings = c('', ' ')) dephosp <- dephosp[,-3] upexpFile <-
  ### system.file('SignedPC', 'dephosphorylates.txt', package='zeptosensPkg') upexp <-
  ### read.csv(upexpFile, sep='\t', header=TRUE, na.strings = c('', ' ')) upexp <- upexp[,-3]
  ### dwnexpFile <- system.file('SignedPC', 'downregulates-expression.txt', package='zeptosensPkg')
  ### dwnexp <- read.csv(dwnexpFile, sep='\t', header=TRUE, na.strings = c('', ' ')) dwnexp <-
  ### dwnexp[,-3]

  #    results <- downloadSignedPC(forceCache=TRUE)
  #    results <- read.table(file="results_network1.txt",header=T)
  results <- read.table(system.file("extdata", "filteredSignedPc.txt", 
                                    package = "zeptosensPkg"), 
                        sep = "\t", 
                        header = TRUE, fill = TRUE)
  #    write.table(results, file="results_network1.txt",quote=F)
  dephosp <- paxtoolsr:::filterSif(results, interactionTypes = "dephosphorylates")
  phosp <- paxtoolsr:::filterSif(results, interactionTypes = "phosphorylates")
  dwnexp <- paxtoolsr:::filterSif(results, interactionTypes = "downregulates-expression")
  upexp <- paxtoolsr:::filterSif(results, interactionTypes = "upregulates-expression")

  # NOTE: SIF has interaction type as column 2, edgelists (like distances) do
  # not have this, so convert the SIF to an edgelist
  dephosp <- dephosp[, c(1, 3)]
  phosp <- phosp[, c(1, 3)]
  dwnexp <- dwnexp[, c(1, 3)]
  upexp <- upexp[, c(1, 3)]

  # only concentration nodes are included in up & downregulation
  # mabToGenes_c <- mabToGenes[which(mabToGenes$Effect == "c"), ]

  # define wk
  wk <- matrix(0, ncol = nProt, nrow = nProt, 
               dimnames = list(colnames(proteomicResponses), colnames(proteomicResponses))) 
  wks <- matrix(0, ncol = nProt, nrow = nProt, 
                dimnames = list(colnames(proteomicResponses), colnames(proteomicResponses))) 

  # upregulation expression, wk=1
  # upexpGene1 <- pmatch(upexp[, 1], mabToGenes_c[measured_genes, 4], duplicates.ok = TRUE)
  # upexpGene2 <- pmatch(upexp[, 3], mabToGenes_c[measured_genes, 4], duplicates.ok = TRUE)
  # upexpGene <- cbind(upexpGene1, upexpGene2)

  # Define genes by their effects
  tmpIdxC <- intersect(idxAbMap, which(mabToGenes$Effect == "c"))
  tmpIdxAC <- intersect(idxAbMap, which(mabToGenes$Effect != "i"))
  tmpIdxAI <- intersect(idxAbMap, which(mabToGenes$Effect != "c"))
  tmpIdxA <- intersect(idxAbMap, which(mabToGenes$Effect == "a"))

  tmpGenesC <- mabToGenes[tmpIdxC, 4]
  tmpGenesA <- mabToGenes[tmpIdxAC, 4]
  tmpGenesD <- mabToGenes[tmpIdxAI, 4]
  tmpGenesAo <- mabToGenes[tmpIdxA, 4]

  names(tmpGenesC) <- mabToGenes[tmpIdxC, 1]
  names(tmpGenesA) <- mabToGenes[tmpIdxAC, 1]
  names(tmpGenesD) <- mabToGenes[tmpIdxAI, 1]
  names(tmpGenesAo) <- mabToGenes[tmpIdxA, 1]

  # only concentration and act. phospho nodes are included in up & downregulation
  upexpGene <- matchGenesToEdgelist(
    genes1 = tmpGenesA, genes2 = tmpGenesC, annotEdgelist = upexp,
    antibodyVec = colnames(proteomicResponses), useAnnot = FALSE, verbose = verbose
  )

  for (i in 1:length(upexp[, 1])) {
    wk[upexpGene[i, 1], upexpGene[i, 2]] <- 1
    wks[upexpGene[i, 1], upexpGene[i, 2]] <- 1

    #        print(upexpGene[i, 1])
  }

  # downregulation expression, wk=-1
  dwnexpGene <- matchGenesToEdgelist(
    genes1 = tmpGenesA, genes2 = tmpGenesC, annotEdgelist = dwnexp,
    antibodyVec = colnames(proteomicResponses), useAnnot = FALSE, verbose = verbose
  )
  # cov318 results in 15

  for (i in 1:length(dwnexp[, 1])) {
    wk[dwnexpGene[i, 1], dwnexpGene[i, 2]] <- -1
    wks[dwnexpGene[i, 1], dwnexpGene[i, 2]] <- -1
  }

  # phosphorylates wk=1 only active and concentration states are upstream
  # mabToGenes_a <- mabToGenes[which(mabToGenes$Effect != "i"), ]
  # mabToGenes_d <- mabToGenes[which(mabToGenes$Sites != "c"), ]
  # phosGene1 <- pmatch(phosp[, 1], mabToGenes_a[measured_genes, 4], duplicates.ok = TRUE)
  # phosGene2 <- pmatch(phosp[, 3], mabToGenes_d[measured_genes, 4], duplicates.ok = TRUE)
  # phosGene <- cbind(phosGene1, phosGene2)

  phosGene <- matchGenesToEdgelist(
    genes1 = tmpGenesAo, genes2 = tmpGenesD, annotEdgelist = phosp,
    antibodyVec = colnames(proteomicResponses), useAnnot = FALSE, verbose = verbose
  )
  # cov318 13 results

  for (i in 1:length(phosGene[, 1])) {
    wk[phosGene[i, 1], phosGene[i, 2]] <- 1
    wks[phosGene[i, 1], phosGene[i, 2]] <- 2
  }

  # dephosphorylates wk=-1 only active and concentration states are upstream
  # mabToGenes_a <- mabToGenes[which(mabToGenes$Effect != "i"), ]
  # mabToGenes_d <- mabToGenes[which(mabToGenes$Sites != "c"), ]
  # dephosGene1 <- pmatch(dephosp[, 1], mabToGenes_a[measured_genes, 4], duplicates.ok = TRUE)
  # dephosGene2 <- pmatch(dephosp[, 3], mabToGenes_d[measured_genes, 4], duplicates.ok = TRUE)
  # dephosGene <- cbind(dephosGene1, dephosGene2)

  dephosGene <- matchGenesToEdgelist(
    genes1 = tmpGenesA, genes2 = tmpGenesD, annotEdgelist = dephosp,
    antibodyVec = colnames(proteomicResponses), useAnnot = FALSE, verbose = verbose
  )

  for (i in 1:length(dephosGene[, 1])) {
    wk[dephosGene[i, 1], dephosGene[i, 2]] <- -1
    wks[dephosGene[i, 1], dephosGene[i, 2]] <- -2
  }
  inter <- (which(wk != 0, arr.ind = T))
  print(inter)
  networks <- list(wk = wk, wks = wks, distInd = distInd, inter = inter)

  return(networks)
}
