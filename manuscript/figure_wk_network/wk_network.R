library(igraph)

mat <- read.table("~/Downloads/wk_x.txt", sep="\t", header=TRUE)

tmpEdgelist <- NULL

for(i in 1:nrow(mat)) {
    #i <- 9
    
    idx <- which(mat[i,] != 0)
    t1 <- colnames(mat)[idx]
    
    if(length(t1) > 0) {
        for(j in 1:length(t1)) {
            t2 <- data.frame(PARTICIPANT_A=rownames(mat)[i], PARTICIPANT_B=t1[j], stringsAsFactors = FALSE)
            tmpEdgelist <- rbind(tmpEdgelist, t2)  
        }
    }
}

i <- 1 

antibodyMapFile <- system.file("targetScoreData", "antibodyMap.txt", package = "zeptosensPkg")      
antibodyMap <- read.table(antibodyMapFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

edgelist <- tmpEdgelist

for(i in 1:nrow(edgelist)) {
    for(j in 1:ncol(edgelist)) {
        cat("IJ: ", i, "::", j, "\n")
        
        idx <- which(grepl(tmpEdgelist[i, j], antibodyMap$AntibodyLabel))
        genes <- antibodyMap$Gene_Symbol[idx]
        
        for(k in 1:length(genes)) {
            edgelist[i, j] <- genes[k]         
        }
    }
}

e2 <- cbind(edgelist, int="hey")
e2 <- unique(e2)
write.table(e2, file="~/Downloads/tmp.txt", sep="\t", row.names=FALSE, quote=FALSE)

library(paxtoolsr)

pcFile <- file.path(.lmp, "gene_set_pathway_analysis", "data", "pc8.rds")

if(file.exists(pcFile)) {
    sif <- readRDS(pcFile)
} else {
    sif <- downloadPc2("PathwayCommons.8.All.EXTENDED_BINARY_SIF.hgnc.txt.gz", version="8")
    saveRDS(sif, pcFile)
}

a <- filterSif(sif$edges, edgelist=edgelist, edgelistCheckReverse=FALSE)
b1 <- downloadSignedPC()
#b2 <- filterSif(b1, edgelist = edgelist)

library(data.table)
b <- a[, c(PARTICIPANT_A, PARTICIPANT_B)]
b2 <- a[, c(1,3), with=FALSE]
setDF(b2)
b3 <- unique(b2)
length(unique(c(b3[,1], b3[,2])))

       
a1 <- downloadSignedPC(forceCache = TRUE)
a2 <- filterSif(a1, ids="MAPK1")





