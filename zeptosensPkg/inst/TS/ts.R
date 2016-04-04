library(zeptosensPkg)
library(zeptosensUtils)
library(paxtoolsr)

nDose=10
nProt=152
maxDist=1 # changing this value requires additional work to compute product(wk). This is not a priority
cellLine="cov318"

#read proteomic response
inputFile <- file.path("inst", "targetScoreData", paste0(cellLine,"_jq1_ave.txt"))
x <- read.table(inputFile, header=TRUE)
#x(dose,prot)
rownames(x) <- x[,1]
x <-x[,-1:-2]
proteomicResponses <- x

targetScoreOutputFile <-"tso.txt"
matrixWkOutputFile <- "wk.txt"
nPerm=20

#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose, nProt, proteomicResponses, 
                          maxDist=1, nPerm,cellLine, 
                          targetScoreOutputFile="~/zeptosenspkg/zeptosenspkg/inst/TS/tso.txt", 
                          matrixWkOutputFile="~/zeptosenspkg/zeptosenspkg/inst/TS/wk.txt",
                          targetScoreQValueFile="~/zeptosenspkg/zeptosenspkg/inst/TS/q.txt", 
                          targetScoreDoseFile="~/zeptosenspkg/zeptosenspkg/inst/TS/tsd.txt",
                          verbose=TRUE)

# 1200 interactions
idx <- which(dist[,3] == 1)
dist[idx, ]

genes <- mab_to_genes[measured_genes, 4]

# 114 genes represented
g2 <- unique(c(dist[idx,1], dist[idx,2]))

# 106 genes from 152 genes represented
length(which(!is.na(pmatch(g2, genes))))
pmatch(g2, genes)

# Interactions without NAs (involve genes of interest)
i1 <- which(complete.cases(dist_list))
# Right dist
i2 <- which(dist_list == 1)
intersect(i1, i2)
[1] 9401 9402 9403 9404 9405 9406 9407 9408 9409 9410 9411 9412 9413 9414 9415 9416 9417 9419 9421 9422 9423 9424
[23] 9425 9426 9427 9428 9429 9430 9431 9432 9433 9434 9435 9436 9437 9438 9439 9440 9441 9442 9443 9444 9445 9446
[45] 9447 9448 9449 9450 9451 9452 9453 9454 9455 9456 9457 9458 9459 9460 9461 9462 9463 9464 9465 9466 9467 9468
[67] 9471 9473 9474 9475 9476 9477 9478 9479 9480 9481 9482 9483 9484 9485 9486 9487 9488 9491 9492 9493 9494 9495
[89] 9496 9497 9498 9499 9500 9501 9502 9503 9504 9505 9506 9507 9508 9509
Browse[2]> length(i1)
[1] 9205
Browse[2]> length(i2)
[1] 1400
Browse[2]> length(i2) / nrow(dist_list)
[1] 0.1320879
Browse[2]> length(intersect(i1, i2)) / nrow(dist_list)
[1] 0.009623549

# 120 to 130 genes represented
which(mab_to_genes[measured_genes, 4] %in% dist[, 1])
which(mab_to_genes[measured_genes, 4] %in% dist[, 2])


