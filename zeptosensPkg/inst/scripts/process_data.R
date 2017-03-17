options("digits"=3)
nrep=2
ncond=3
ntime=4
nab=72
ctr=ncond*ntime*nrep
ct=ncond*ntime
data <- read.table("a2058.txt")
annot <- read.table("ab_index.txt")
normal <- read.table("cellmix.txt")
for (i in 1:length(data[,1])) {
	if(data[i,5] == 'Poor') { 
data[i,3] = 0
data[i,4] = 0
}
}
data.matrix <- as.matrix(data)
annot.matrix <- as.matrix(annot)
datar <- array(dim=c(nab*ncond*ntime,7))

for(ab in 1:nab){
	for (k in 1:ct) {
#		print(k)
	datar[ct*(ab-1)+k,1]=data.matrix[ctr*(ab-1)+k,2]
	datar[ct*(ab-1)+k,2]=annot.matrix[ab,1]
	datar[ct*(ab-1)+k,3]=annot.matrix[ab,2]	
	datar[ct*(ab-1)+k,4]=annot.matrix[ab,3]	
	datar[ct*(ab-1)+k,5]=annot.matrix[ab,4]		
#	datar[ct*(ab-1)+k,5]=data.matrix[ctr*(ab-1)+k,2]	
#	datar[ct*(ab-1)+k,6]=data.matrix[ctr*(ab-1)+k,3]
#	datar[ct*(ab-1)+k,7]=data.matrix[ctr*(ab-1)+ct+k,3]
	datar[ct*(ab-1)+k,6]=data[ctr*(ab-1)+k,4]#/(normal[ab,4])
	datar[ct*(ab-1)+k,7]=data[ctr*(ab-1)+ct+k,4]#/(normal[ab,4])
	}

}
#print(datar[,4])
write.table(datar, file="a2058_replicates_nn.txt", col.names=F, row.names=F, sep = "	",quote=FALSE)
g_datar=datar[which(as.numeric(datar[,6]) > 0 & as.numeric(datar[,7]) > 0 ),]
#print(g_datar)
cc=cor(as.numeric(g_datar[,6]),as.numeric(g_datar[,7]))

pdf("replicates.pdf")
plot(g_datar[,6],g_datar[,7],xlab="replicate 1",ylab="replicate 2")
title(main=paste("a2058 replicates \nCC=",cc))
dev.off()
pdf("replicates_2.pdf")
plot(g_datar[,6],g_datar[,7],xlab="replicate 1",ylab="replicate 2",xlim=c(0,2),ylim=c(0,2))
title(main=paste("a2058 replicates\nCC=",cc))
dev.off()
rm(list=ls())

