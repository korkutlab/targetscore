#correlations w/ prot vs cv
library(Hmisc)
cov318_cv <- read.table(file="cov318_cv.txt",header=TRUE)
cov318_pr <- read.table(file="cov318_jq1_ave.txt",header=TRUE)
print(cov318_cv)
file.remove("cov318_spR.txt")
#print(colnames(cov318_pr[3]))
for(k in 3:length(cov318_pr) ){
	t1 <-matrix(c(as.numeric(cov318_pr[,k]),as.numeric(cov318_cv[,3])), nrow=10,ncol=2)	
#	print(t1)
tc <- rcorr(t1, type="spearman")
print(c(colnames(cov318_pr[k]),tc$r[2,1],tc$P[2,1]))
tp <-matrix(c(colnames(cov318_pr[k]),tc$r[2,1],tc$P[2,1]),nrow=1,ncol=3)
write.table(tp,file="cov318_spR.txt",append=T,row.names=F,col.names=F,quote=F)
#write.table()
#print(c(colnames(cov318_pr[k]),tc[1,2]))
}
nox <- read.table(file="cov318_spR.txt",header=FALSE)
#x2 <- as.numeric(nox[,2])
#print(nox[,2])
sorti <-nox[order(nox[,2],decreasing=T),]
#incl_list()
write.table(sorti,file="cov318_spr_sorted.txt",row.names=F,col.names=F,quote=F)
#rankedcor <- nox
print(sorti)
#print(rankedcor)
rm(list=ls()) 
