library(readxl)

dat <- read_excel("./process_mdacc_protein_list/protein_lists.xlsx", sheet=1, col_names=FALSE)

abIds <- NULL
simplifyFlag <- FALSE

for (i in 1:nrow(dat)) {
  for(j in 1:ncol(dat)) {
    if(!is.na(dat[i,j])) {
      
      tmp <- dat[i,j]
      
      if(simplifyFlag) {
        # Remove suffixes
        tmp <- gsub('\\(.*\\)', '', tmp)
        tmp <- gsub('\\.txt', '', tmp)
        tmp <- gsub('_GBL.*$', '', tmp) 
        tmp <- gsub('_GBL.*$', '', tmp) 
        tmp <- gsub('_R[Ss].*$', '', tmp) 
        
        tmp <- gsub('[\\.-]w?[RGM][\\.-][QVEC].*$', '', tmp) 
        tmp <- gsub('\\.R$', '', tmp) 
        
        # To uppercase
        tmp <- toupper(tmp)
        
        # Change dashes
        tmp <- gsub('\\.', '_', tmp) 
        tmp <- gsub('-', '_', tmp) 
        
        # For comparison
        tmp <- gsub('_', '', tmp) 
        tmp <- gsub('_', '', tmp) 
      }

      abIds <- c(abIds, tmp)
    }
  }
}

tmp <- sort(unique(abIds))

df <- data.frame("AntibodyLabel"=tmp,
                 "Source"=rep("MDACC", length(tmp)),
                 "NodeName"=rep("", length(tmp)),
                 "Gene_Symbol"=rep("", length(tmp)),
                 "Sites"=rep("", length(tmp)),
                 "Effect"=rep("", length(tmp)))


if(simplifyFlag) {
  write.table(df, "./process_mdacc_protein_list/unique_abs.txt", col.names = FALSE, row.names = FALSE, quote=FALSE)  
} else {
  write.table(df, "./process_mdacc_protein_list/unique_abs_full.txt", col.names = FALSE, row.names = FALSE, quote=FALSE)
}


