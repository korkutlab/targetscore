#correlations w/ prot vs cv
library(Hmisc)

rm(list = ls())

basePath <- "figure_spearman_heatmaps/spearmanData/"
cellLines <- c("cov318", "igrov1", "ovcar3", "ovcar4")

for(cellLine in cellLines) {
  #DEBUG 
  #cellLine <- "cov318"
  
  cellViabilityValues <- read.table(file = paste0(basePath, cellLine, "_cv.txt"), header = TRUE)
  protValues <- read.table(file = paste0(basePath, cellLine, "_jq1_ave.txt"), header = TRUE)
  
  # DEBUG
  #print(cellViabilityValues)
  
  file.remove(paste0(basePath, cellLine, "_spR.txt"))
  #print(colnames(protValues[3]))
  
  for (k in 3:length(protValues)) {
    #DEBUG 
    #k <- 3
    t1 <-
      matrix(c(as.numeric(protValues[, k]), as.numeric(cellViabilityValues[, 3])),
             nrow = 10,
             ncol = 2)
    #	print(t1)
    tc <- rcorr(t1, type = "spearman")
    #print(c(colnames(protValues[k]), tc$r[2, 1], tc$P[2, 1]))
    tp <- matrix(c(colnames(protValues[k]), tc$r[2, 1], tc$P[2, 1]), nrow = 1, ncol = 3)
    write.table(
      tp,
      file = paste0(basePath, cellLine, "_spR.txt"),
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    #write.table()
    #print(c(colnames(protValues[k]),tc[1,2]))
  }
  
  nox <- read.table(file = paste0(basePath, cellLine, "_spR.txt"), header = FALSE)
  colnames(nox) <- c("abName", "corr", "pval")
  #x2 <- as.numeric(nox[,2])

  # ADJUST PVALUES
  nox[, "adjPval"] <- p.adjust(nox$pval, method="fdr")

  # REORDER VALUES
  #print(nox[,2])
  sorti <- nox[order(nox[, "corr"], decreasing = TRUE), ]
  #incl_list()
  write.table(
    sorti,
    file = paste0(basePath, cellLine, "_spr_sorted.txt"),
    sep="\t",
    row.names = FALSE,
    quote = FALSE
  )
  #rankedcor <- nox
  #print(sorti)
  #print(rankedcor)
}
