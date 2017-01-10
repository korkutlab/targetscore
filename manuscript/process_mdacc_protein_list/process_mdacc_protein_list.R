library(readxl)

dat <- read_excel("./process_mdacc_protein_list/protein_lists.xlsx", sheet=1, col_names=FALSE)

abIds <- NULL

for (i in 1:nrow(dat)) {
  for(j in 1:ncol(dat)) {
    if(!is.na(dat[i,j])) {
      
      tmp <- dat[i,j]
      
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
      
      abIds <- c(abIds, tmp)
    }
  }
}

tmp <- sort(unique(abIds))

write.table(tmp, "./process_mdacc_protein_list/unique_abs.txt", col.names = FALSE, row.names = FALSE, quote=FALSE)

