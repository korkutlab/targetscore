#' Target Score Pilot Run
#' 
#' @param nDose TBA
#' @param nProt TBA
#' @param proteomicResponses TBA
#' @param maxDist TBA (default: 1)
#' @param cellLine TBA 
#' @param targetScoreOutputFile a filename to write target score results (default: NULL)
#' @param matrixWkOutputFile TBA 
#' 
#' @details 
#' data: multiple dose single drug perturbation
#' ts: integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))
#' missing: For phosp and dephosp based wk, there is no "exact match" between known and measured phospho-sites
#' 
#' @examples 
#' 
#' @concept zeptosensPkg
#' @export
getTargetScore <- function(nDose, nProt, proteomicResponses, maxDist=1,  
                           cellLine, targetScoreOutputFile=NULL, matrixWkOutputFile=NULL) {
    # LOAD INTERNAL DATA ----
    #read function score
    fsFile <- system.file("targetScoreData", "fs.txt", package="zeptosensPkg")
    fs <- read.table(fsFile, header=TRUE)
    #print(fs)
    
    #match Ab names to gene names & posttranslational modifications
    antibodyMapFile <- system.file("targetScoreData", "antibodyMap.txt", package="zeptosensPkg")
    mab_to_genes <- read.table(antibodyMapFile, sep="\t", header=TRUE)
    #mab_to_genes
    
    #pathway distance matrix
    distFile <- system.file("targetScoreData", "antibodyMap.txt", package="zeptosensPkg")
    dist <- read.table(distFile, sep="\t", header=TRUE)
    #dist
    
    measured_genes <- pmatch(colnames(proteomicResponses),mab_to_genes[,1],duplicates.ok = TRUE)
    mab_to_genes[measured_genes,]
    dist_gene1 <- pmatch(dist[,1],mab_to_genes[measured_genes,4],duplicates.ok = TRUE)
    dist_gene2 <- pmatch(dist[,2],mab_to_genes[measured_genes,4],duplicates.ok = TRUE)
    
    #distance framework
    dist_list <- cbind(dist_gene1,dist_gene2,dist[,3])
    dist_list[is.na(dist_list[,1]),1]<-100
    dist_list[is.na(dist_list[,2]),2]<-100
    dist_list[,1]
    dist_ind <- matrix(Inf,ncol=nProt,nrow=nProt) #dist_ind(upstream,downstream)
    
    for(i in 1:length(dist_list[,1])){
        
        dist_ind[dist_list[i,1],dist_list[i,2]] <- dist[i,3] 
        
        if(dist_ind[dist_list[i,1],dist_list[i,2]] > maxDist){
            dist_ind[dist_list[i,1],dist_list[i,2]] <- Inf
        }
        if(dist_ind[dist_list[i,1],dist_list[i,2]]==0){
            dist_ind[dist_list[i,1],dist_list[i,2]] <- Inf
        }
    }
    
    ###get the network product###
    
    phospFile <- system.file("SignedPC", "phosphorylates.txt", package="zeptosensPkg")
    phosp <- read.csv(phospFile, sep="\t", header=TRUE, na.strings = c("", " "))
    phosp <- phosp[,-3]
    
    dephospFile <- system.file("SignedPC", "dephosphorylates.txt", package="zeptosensPkg")
    dephosp <- read.csv(dephospFile, sep="\t", header=TRUE, na.strings = c("", " "))
    dephosp <- dephosp[,-3]
    
    upexpFile <- system.file("SignedPC", "dephosphorylates.txt", package="zeptosensPkg")
    upexp <- read.csv(upexpFile, sep="\t", header=TRUE, na.strings = c("", " "))
    upexp <- upexp[,-3]
    
    dwnexpFile <- system.file("SignedPC", "downregulates-expression.txt", package="zeptosensPkg")
    dwnexp <- read.csv(dwnexpFile, sep="\t", header=TRUE, na.strings = c("", " "))
    dwnexp <- dwnexp[,-3]
    
    #only concentration nodes are included in up & downregulation
    mab_to_genes_c <- mab_to_genes[which(mab_to_genes$Effect=='c'),]
    
    #define wk
    wk <- matrix(0,ncol=nProt,nrow=nProt) #wk(upstr,downstr)
    #upregulation expression, wk=1
    upexp_gene1 <- pmatch(upexp[,1],mab_to_genes_c[measured_genes,4], duplicates.ok = TRUE)
    upexp_gene2 <- pmatch(upexp[,2],mab_to_genes_c[measured_genes,4], duplicates.ok = TRUE)
    upexp_gene <- cbind(upexp_gene1,upexp_gene2)
    
    for(i in 1:length(upexp[,1])){
        wk[upexp_gene[i,1],upexp_gene[i,2]] = 1
    }
    
    #downregulation expression, wk=-1
    dwnexp_gene1 <- pmatch(dwnexp[,1],mab_to_genes_c[measured_genes,4],duplicates.ok = TRUE)
    dwnexp_gene2 <- pmatch(dwnexp[,2],mab_to_genes_c[measured_genes,4],duplicates.ok = TRUE)
    dwnexp_gene <- cbind(dwnexp_gene1,dwnexp_gene2)
    
    for(i in 1:length(dwnexp[,1])){
        wk[dwnexp_gene[i,1],dwnexp_gene[i,2]] = -1
    }
    
    #phosphorylates wk=1
    #only active and concentration states are upstream
    mab_to_genes_a <- mab_to_genes[which(mab_to_genes$Effect != 'i'),]
    mab_to_genes_d <- mab_to_genes[which(mab_to_genes$Sites != 'c'),]
    phos_gene1 <- pmatch(phosp[,1],mab_to_genes_a[measured_genes,4],duplicates.ok = TRUE) 
    phos_gene2 <- pmatch(phosp[,2],mab_to_genes_d[measured_genes,4],duplicates.ok = TRUE) 
    phos_gene <- cbind(phos_gene1,phos_gene2)
    
    for(i in 1:length(phos_gene[,1])){
        wk[phos_gene[i,1],phos_gene[i,2]] = 1 
    }
    
    #dephosphorylates wk=-1
    #only active and concentration states are upstream
    mab_to_genes_a <- mab_to_genes[which(mab_to_genes$Effect != 'i'),]
    mab_to_genes_d <- mab_to_genes[which(mab_to_genes$Sites != 'c'),]
    dephos_gene1 <- pmatch(dephosp[,1],mab_to_genes_a[measured_genes,4],duplicates.ok = TRUE) 
    dephos_gene2 <- pmatch(dephosp[,2],mab_to_genes_d[measured_genes,4],duplicates.ok = TRUE) 
    dephos_gene <- cbind(dephos_gene1,dephos_gene2)
    
    for(i in 1:length(dephos_gene[,1])){
        wk[dephos_gene[i,1],dephos_gene[i,2]] = -1
    }
    
    if(!is.null(matrixWkOutputFile)) {
        write.table(wk, file=matrixWkOutputFile)
    }
    
    #calculate TS for each dose
    ts <- matrix(0,ncol=nProt,nrow=nDose)
    tsp <- array(0:0,dim=c(nDose,nProt,nProt))
    for(i in 1:nDose) {
        #downstream (target)
        for (j in 1:nProt)  {
            #upstream
            for (k in 1:nProt){
                
                tsp[i,k,j]=(2^-(dist_ind[k,j]))*proteomicResponses[i,k]*wk[k,j]
                
            }
            ts[i,j]=fs[i,2]*(proteomicResponses[i,j]+sum(tsp[i,1:nProt,j]))
        }
    }
    colnames(ts) <- colnames(proteomicResponses)
    rownames(ts) <- rownames(proteomicResponses)
    
    if(!is.null(targetScoreOutputFile)) {
        write.table(ts, file=targetScoreOutputFile)
    }

    results <- list(ts=ts, wk=wk)
    return(results)
}
  