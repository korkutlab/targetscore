library(readxl)

dat <- read_excel("./process_mdacc_protein_list/protein_lists.xlsx", sheet=1, col_names=FALSE)

abIds <- NULL
simpleNames <- NULL 
tmpDf <- NULL

for (i in 1:nrow(dat)) {
  for(j in 1:ncol(dat)) {
    if(!is.na(dat[i,j])) {
      
      tmp <- dat[i,j]
      abIds <- c(abIds, tmp)
      
      # Remove suffixes
      tmp <- gsub('\\(.*\\)', '', tmp)
      tmp <- gsub('\\.txt', '', tmp)
      tmp <- gsub('_GBL.*$', '', tmp) 
      tmp <- gsub('_GBL.*$', '', tmp) 
      tmp <- gsub('_R[Ss].*$', '', tmp) 
      
      tmp <- gsub('[\\.-]w?[RGM][\\.-][QVEC].*$', '', tmp) 
      tmp <- gsub('\\.R$', '', tmp) 
      
      # To uppercase
      #tmp <- toupper(tmp)
      
      # Change dashes
      tmp <- gsub('\\.', '_', tmp) 
      tmp <- gsub('-', '_', tmp) 
      
      # For table
      nodeName <- gsub('_', '', tmp) 
      symbol <- sub("^(.*)_p.*$", "\\1", tmp)
      symbol <- gsub('_', '', symbol) 
      symbol <- toupper(symbol)
      
      sites <- ''
      if(grepl("p[TYS].*$", tmp)) {
        sites <- sub(".*(p[TYS].*)$", "\\1", tmp)
      }

      tmpDf <- rbind(tmpDf, data.frame(AntibodyLabel=dat[i,j], 
                                       Source='MDACC',
                                       NodeName=nodeName, 
                                       Gene_Symbol=symbol, 
                                       Sites=sites, 
                                       Effect='',
                                       FunctionalScore='',
                                       stringsAsFactors=FALSE)) 
    }
  }
}

tmpDf <- unique(tmpDf)
tmpDf <- tmpDf[with(tmpDf, order(AntibodyLabel)), ]

#write.table(tmp, "./process_mdacc_protein_list/unique_abs.txt", col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(tmpDf, "./process_mdacc_protein_list/unique_abs_full.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)



