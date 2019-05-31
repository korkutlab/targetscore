#' Extracted functional score value from COMIC/ONCODRIVE Database. Can be override with Manually set functional score.
#' 
#' @param nProt number of proteins tested within the data.
#' @param proteomicResponses Proteomic responses as drug pertubation.
#' @param antibodyMapFile a listing of antibodies, their associated genes, and modification sites
#' @param fsValueFile a listing of functional scores for each gene Manually set uo for overriding Cosmic Database given value, the modification path. (.txt)
#' @param verbose Defult as Fault.If given True, will print out the gene seq mapped with Antibody Map File.
#' @examples 
#' 
#' @concept zeptosensPkg
#' @export
getFsVals <- function(nProt,proteomicResponses,antibodyMapFile=NULL,fsValueFile=NULL,verbose=F){
    
    # match Ab names to gene names & posttranslational modifications
    if(is.null(antibodyMapFile)) {
        antibodyMapFile <- system.file("targetScoreData", "antibodyMap.txt", package = "zeptosensPkg")      
    }
    mab_to_genes <- read.table(antibodyMapFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    if(verbose) {
        print(mab_to_genes)   
    }
    
    if(nProt != ncol(proteomicResponses)) {
        stop("ERROR: nProt is not equal to proteomicResponses column number")
    }
    
    #Match the protein names in the proteomicresponce with the AntibodyMapfile
    idxAbMap <- which(mab_to_genes[, 1] %in% colnames(proteomicResponses))
    if(length(idxAbMap) < nProt) {
        stop("ERROR: Not all columns in data were matched in antibody map")
    }
    if(length(unique(mab_to_genes[idxAbMap, 1])) != nProt) {
        print(unique(mab_to_genes[idxAbMap, 1]))
        stop("ERROR: Mismatch in the number of selected antibodies and the number of proteomic responses")
    }
    antibodyMapSubset <- mab_to_genes[idxAbMap, ]
    
    mabGenes <- mab_to_genes[idxAbMap, 4]   
    names(mabGenes) <- mab_to_genes[idxAbMap, 1]
    
    #Get FS value
    mabValue <- mab_to_genes[idxAbMap, 6]
    mabFS <- ifelse(mabValue=="a",1,ifelse(mabValue=="i",-1,ifelse(mabValue=="c",1,0))) 
    
    CancerRole <- read.table(system.file("extdata", "Cosmic.txt", package="zeptosensPkg"), sep="\t", header=TRUE, fill=TRUE)    
    Cos_FS<-CancerRole[mabGenes,]$fs
    
    fs_value=mabFS*Cos_FS
    
    fs=data.frame(prot=mab_to_genes[idxAbMap, 1],fs=fs_value)
        
    #Uniqueness of fs value
    prot=as.character(unlist(unique(fs$prot)))
    fs=unique(fs)
        
    #match with antibody seq
    index=match(colnames(proteomicResponses),fs$prot)
    fs=fs[index,]
    
    #Override with Self setting/external fs value
    if(!is.null(fs.override)){
    fs.override=read.table(fsValueFile,sep="\t", header=TRUE, fill=TRUE)
    index=which(fs$prot%in%fs.override$prot)
    fs[index,2]<-fs.override$fs
    }
    
    return(fs)
}