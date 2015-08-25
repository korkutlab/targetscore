library(xlsx)
colClasses <- c("character","character","numeric")
        
tpl_1 <- read.xlsx2("inst/dataInst/R022_A2058_sample_conc_EZQ.xlsx",colClasses=colClasses, stringsAsFactors=F, sheetIndex=1)
tpl_2 <- read.xlsx2("inst/dataInst/R023_Mel133_sample_conc_EZQ.xlsx",colClasses=colClasses, stringsAsFactors=F, sheetIndex=1)
print(tpl_1)

tpl_c <- tpl_1[,3]/tpl_2[,3]
print(tpl_c)
tpl_n <- cbind(tpl_1[,1],tpl_1[,2],tpl_c)
print(tpl_n)
colnames(tpl_n) <- c("sample_ID","sample_name","tpl_norm_const")

write.table(tpl_n,file="inst/dataInst/tpl_norm.txt",quote=F,row.names=F, sep="\t")
rm(list=ls())
