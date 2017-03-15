#network from signedPC
#' @param nProt TBA
#' @param proteomicResponses TBA
#' @param maxDist TBA (default: 1)
#' @param antibodyMapFile a listing of antibodies, their associated genes, and modification sites
#' @param distFile A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' 
#' 
#' @examples 
#' 
#' @concept zeptosensPkg
#' @export
network <- function(nProt,proteomicResponses,maxDist,antibodyMapFile=NULL,distFile=NULL,verbose=F){

    
    # match Ab names to gene names & posttranslational modifications
    if(is.null(antibodyMapFile)) {
        antibodyMapFile <- system.file("targetScoreData", "antibodyMap.txt", package = "zeptosensPkg")      
    }
    mab_to_genes <- read.table(antibodyMapFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    if(verbose) {
        print(mab_to_genes)   
    }
    
    # pathway distance matrix
    if(is.null(distFile)) {
        distFile <- system.file("targetScoreData", "distances.txt", package = "zeptosensPkg")   
    }
    
    tmpDist <- read.table(distFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    if(nProt != ncol(proteomicResponses)) {
        stop("ERROR: nProt is not equal to proteomicResponses column number")
    }
    
    # Filter dist to only keep those with a distance less than maxDist
    idx <- which(tmpDist[,3] <= maxDist)
    dist <- tmpDist[idx,]
    # dist
    idxAbMap <- which(mab_to_genes[, 1] %in% colnames(proteomicResponses))
    #    print(mab_to_genes[, 1])
    #    print(unique(mab_to_genes[idxAbMap,]))
    #    if(length(idxAbMap) < nProt) {
    #    print(length(unique(mab_to_genes[idxAbMap,1])))
    if(length(idxAbMap) < nProt) {
        stop("ERROR: Not all columns in data were matched in antibody map")
    }
    #    print((unique(mab_to_genes[idxAbMap, 1])))
    if(length(unique(mab_to_genes[idxAbMap, 1])) != nProt) {
        stop("ERROR: Mismatch in the number of selected antibodies and the number of proteomic responses")
    }
    
    # Used by matchGenesToEdgelist to account for the cases where multiple entries for the same antibody
    # exist in the antibody map
    antibodyMapSubset <- mab_to_genes[idxAbMap, ]
    
    mabGenes <- mab_to_genes[idxAbMap, 4]   
    names(mabGenes) <- mab_to_genes[idxAbMap, 1]
    
    dist_list <- matchGenesToEdgelist(genes1 = mabGenes, genes2 = NULL, annotEdgelist = dist, 
                                      antibodyVec=colnames(proteomicResponses), useAnnot=TRUE, verbose = TRUE)
    
    
    # # Get interactions for measured genes
    # dist_gene1 <- pmatch(dist[, 1], mab_to_genes[measured_genes, 4], duplicates.ok = TRUE)
    # dist_gene1Name <- mab_to_genes[measured_genes, 4][dist_gene1]
    # dist_gene2 <- pmatch(dist[, 2], mab_to_genes[measured_genes, 4], duplicates.ok = TRUE)
    # dist_gene2Name <- mab_to_genes[measured_genes, 4][dist_gene2]
    # 
    # # distance framework
    # dist_list <- data.frame(dist_gene1=dist_gene1, dist_gene2=dist_gene2, dist=dist[, 3], 
    #                         dist_gene1Name=dist_gene1Name, dist_gene2Name=dist_gene2Name, stringsAsFactors = FALSE)
    # dist_list[is.na(dist_list[, 1]), 3] <- 100
    # dist_list[is.na(dist_list[, 2]), 3] <- 100
    # dist_list[, 1]
    dist_ind <- matrix(Inf, ncol = nProt, nrow = nProt, dimnames=list(colnames(proteomicResponses), colnames(proteomicResponses)))  #dist_ind(upstream,downstream)
    
    for (i in 1:length(dist_list[, 1])) {
        dist_ind[dist_list[i, 1], dist_list[i, 2]] <- dist_list[i, 3]
        
        if (dist_ind[dist_list[i, 1], dist_list[i, 2]] > maxDist) {
            dist_ind[dist_list[i, 1], dist_list[i, 2]] <- Inf
        }
        
        if (dist_ind[dist_list[i, 1], dist_list[i, 2]] == 0) {
            dist_ind[dist_list[i, 1], dist_list[i, 2]] <- Inf
        }
    }
    write.table(dist_ind,file="dist_ind.txt", quote=F)
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
    results <- read.table(system.file("extdata", "filteredSignedPc.txt", package="zeptosensPkg"), sep="\t", header=TRUE, fill=TRUE)    
    #    write.table(results, file="results_network1.txt",quote=F)
    dephosp <- filterSif(results, interactionTypes="dephosphorylates")
    phosp <- filterSif(results, interactionTypes="phosphorylates")
    dwnexp <- filterSif(results, interactionTypes="downregulates-expression")
    upexp <- filterSif(results, interactionTypes="upregulates-expression")
    
    # NOTE: SIF has interaction type as column 2, edgelists (like distances) do 
    # not have this, so convert the SIF to an edgelist
    dephosp <- dephosp[,c(1,3)]
    phosp <- phosp[,c(1,3)]
    dwnexp <- dwnexp[,c(1,3)]
    upexp <- upexp[,c(1,3)]
    
    # only concentration nodes are included in up & downregulation
    #mab_to_genes_c <- mab_to_genes[which(mab_to_genes$Effect == "c"), ]
    
    # define wk
    wk <- matrix(0, ncol = nProt, nrow = nProt, dimnames = list(colnames(proteomicResponses), colnames(proteomicResponses)))  #wk(upstr,downstr)
    wks <- matrix(0, ncol = nProt, nrow = nProt, dimnames = list(colnames(proteomicResponses), colnames(proteomicResponses)))  #wk(upstr,downstr)
    
    # upregulation expression, wk=1
    # upexp_gene1 <- pmatch(upexp[, 1], mab_to_genes_c[measured_genes, 4], duplicates.ok = TRUE)
    # upexp_gene2 <- pmatch(upexp[, 3], mab_to_genes_c[measured_genes, 4], duplicates.ok = TRUE)
    # upexp_gene <- cbind(upexp_gene1, upexp_gene2)
    
    # Define genes by their effects
    tmpIdxC <- intersect(idxAbMap, which(mab_to_genes$Effect == "c"))
    tmpIdxAC <- intersect(idxAbMap, which(mab_to_genes$Effect != "i"))
    tmpIdxAI <- intersect(idxAbMap, which(mab_to_genes$Effect != "c"))
    tmpIdxA <- intersect(idxAbMap, which(mab_to_genes$Effect == "a"))
    
    tmpGenesC <- mab_to_genes[tmpIdxC, 4]
    tmpGenesA <- mab_to_genes[tmpIdxAC, 4]
    tmpGenesD <- mab_to_genes[tmpIdxAI, 4]
    tmpGenesAo <- mab_to_genes[tmpIdxA, 4]
    
    names(tmpGenesC) <- mab_to_genes[tmpIdxC, 1]
    names(tmpGenesA) <- mab_to_genes[tmpIdxAC, 1]
    names(tmpGenesD) <- mab_to_genes[tmpIdxAI, 1]
    names(tmpGenesAo) <- mab_to_genes[tmpIdxA, 1]
    
    # only concentration and act. phospho nodes are included in up & downregulation
    upexp_gene <- matchGenesToEdgelist(genes1=tmpGenesA, genes2=tmpGenesC, annotEdgelist=upexp, 
                                       antibodyVec=colnames(proteomicResponses), useAnnot=FALSE, verbose=verbose)
    
    for (i in 1:length(upexp[, 1])) {
        wk[upexp_gene[i, 1], upexp_gene[i, 2]] <- 1
        wks[upexp_gene[i, 1], upexp_gene[i, 2]] <- 1
        
        #        print(upexp_gene[i, 1])
    }
    
    # downregulation expression, wk=-1
    dwnexp_gene <- matchGenesToEdgelist(genes1=tmpGenesA, genes2=tmpGenesC, annotEdgelist=dwnexp, 
                                        antibodyVec=colnames(proteomicResponses), useAnnot=FALSE, verbose=verbose)
    # cov318 results in 15
    
    for (i in 1:length(dwnexp[, 1])) {
        wk[dwnexp_gene[i, 1], dwnexp_gene[i, 2]] <- -1
        wks[dwnexp_gene[i, 1], dwnexp_gene[i, 2]] <- -1
    }
    
    # phosphorylates wk=1 only active and concentration states are upstream
    # mab_to_genes_a <- mab_to_genes[which(mab_to_genes$Effect != "i"), ]
    # mab_to_genes_d <- mab_to_genes[which(mab_to_genes$Sites != "c"), ]
    # phos_gene1 <- pmatch(phosp[, 1], mab_to_genes_a[measured_genes, 4], duplicates.ok = TRUE)
    # phos_gene2 <- pmatch(phosp[, 3], mab_to_genes_d[measured_genes, 4], duplicates.ok = TRUE)
    # phos_gene <- cbind(phos_gene1, phos_gene2)
    
    phos_gene <- matchGenesToEdgelist(genes1=tmpGenesAo, genes2=tmpGenesD, annotEdgelist=phosp, 
                                      antibodyVec=colnames(proteomicResponses), useAnnot=FALSE, verbose=verbose)
    #cov318 13 results
    
    for (i in 1:length(phos_gene[, 1])) {
        wk[phos_gene[i, 1], phos_gene[i, 2]] <- 1
        wks[phos_gene[i, 1], phos_gene[i, 2]] <- 2
        
    }
    
    # dephosphorylates wk=-1 only active and concentration states are upstream
    # mab_to_genes_a <- mab_to_genes[which(mab_to_genes$Effect != "i"), ]
    # mab_to_genes_d <- mab_to_genes[which(mab_to_genes$Sites != "c"), ]
    # dephos_gene1 <- pmatch(dephosp[, 1], mab_to_genes_a[measured_genes, 4], duplicates.ok = TRUE)
    # dephos_gene2 <- pmatch(dephosp[, 3], mab_to_genes_d[measured_genes, 4], duplicates.ok = TRUE)
    # dephos_gene <- cbind(dephos_gene1, dephos_gene2)
    
    dephos_gene <- matchGenesToEdgelist(genes1=tmpGenesA, genes2=tmpGenesD, annotEdgelist=dephosp, 
                                        antibodyVec=colnames(proteomicResponses), useAnnot=FALSE, verbose=verbose)

    for (i in 1:length(dephos_gene[, 1])) {
        wk[dephos_gene[i, 1], dephos_gene[i, 2]] <- -1
        wks[dephos_gene[i, 1], dephos_gene[i, 2]] <- -2
        
    }  
    inter <-(which(wk!=0,arr.ind=T))
    print(inter)
    networks  <- list(wk = wk,wks=wks, dist_ind=dist_ind, inter=inter)
    
    return(networks)
}


