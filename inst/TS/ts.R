# target score pilot run
#anil/augustin/ozgun 11/15
#data: multiple dose single drug perturbation
#ts = integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))
#missing: For phosp and dephosp based wk, there is no "exact match" between known and measured phospho-sites
ndose=11
nprot=154
maxdist=1 # changing this value requires additional work to compute product(wk). This is not a priority
data_dir="~/zeptosenspkg/TS_data/"
cell_line="cov318"

#read proteomic response
x <- read.table(paste(data_dir,cell_line,"_jq1_ave.txt",sep=""),header=T)
rownames(x) <- x[,1]
x <-x[,-1:-2]
#read function score
fs <- read.table(paste(data_dir,"fs.txt",sep=""),header=T)
print(fs)
#match Ab names to gene names & posttranslational modifictions
mab_to_genes <- read.table(paste(data_dir,"antibodyMap.txt",sep=""),sep="\t",header=T)
mab_to_genes
#pathway distance matrix
dist <- read.table(paste(data_dir,"distances.txt",sep=""),sep="\t",header=T)
dist
measured_genes <- pmatch(colnames(x),mab_to_genes[,1],duplicates.ok = TRUE)
mab_to_genes[measured_genes,]
dist_gene1 <- pmatch(dist[,1],mab_to_genes[measured_genes,4],duplicates.ok = TRUE)
dist_gene2 <- pmatch(dist[,2],mab_to_genes[measured_genes,4],duplicates.ok = TRUE)

#distance framework
dist_list <- cbind(dist_gene1,dist_gene2,dist[,3])
dist_list[is.na(dist_list[,1]),1]<-Inf
dist_list[is.na(dist_list[,2]),2]<-Inf
dist_list[,1]
dist_ind <- matrix(0,ncol=nprot,nrow=nprot)

for(i in 1:length(dist_list[,1])){

  dist_ind[dist_list[i,1],dist_list[i,2]] <- dist[i,3] 
  }
###get the network product###

phosp <- read.csv(paste(data_dir,"phosphorylates.txt",sep=""),sep="\t",header=T,na.strings = c("", " "))
phosp <- phosp[,-3]
dephosp <- read.csv(paste(data_dir,"dephosphorylates.txt",sep=""),sep="\t",header=T,na.strings = c("", " "))
dephosp <- dephosp[,-3]
upexp <- read.csv(paste(data_dir,"upregulates-expression.txt",sep=""),sep="\t",header=T,na.strings = c("", " "))
upexp <- upexp[,-3]
dwnexp <- read.csv(paste(data_dir,"downregulates-expression.txt",sep=""),sep="\t",header=T,na.strings = c("", " "))
dwnexp <- dwnexp[,-3]
#only concentration nodes are included in up & downregulation
mab_to_genes_c <- mab_to_genes[which(mab_to_genes$Effect=='c'),]
#define wk
wk <- matrix(0,ncol=nprot,nrow=nprot)
#upregulation expression, wk=1
upexp_gene1 <- pmatch(upexp[,1],mab_to_genes_c[measured_genes,4],duplicates.ok = TRUE)
upexp_gene2 <- pmatch(upexp[,2],mab_to_genes_c[measured_genes,4],duplicates.ok = TRUE)
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
write.table(wk,file="wk.txt")
#IN PROGRESS!
#calculate TS for each dose
TS <- matrix(0,ncol=ndose,nrow=nprot)

for(i in 1:ndose){
  for (j in 1:nprot)
    ts[i,j]=fs[i,1]*(x[i,1]+
}

rm(list = ls()) 

#permutation 
  