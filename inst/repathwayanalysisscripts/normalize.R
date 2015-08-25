options("digits"=3)
nrep=2
ncond=4
ntime=3
nab=72
ctr=ncond*ntime*nrep
ct=ncond*ntime
data <- read.table("a2058.txt")
antib <- read.table("ab_index.txt")
print(data[4,2])
exper<-(strsplit(as.character(data[,2]),"_"))

for (i in 1:length(data[,1])) {
	if(data[i,5] == 'Poor') { 
data[i,3] = 0
data[i,4] = 0
}
}
normal <- read.table("cellmix.txt")
#normal_n <- as.matrix(normal)
x <- array(dim=c(nab,ntime,ncond,7))
data_list=array(dim=c(nab*ntime*ncond,7))
print(x[1,1,1,1])

for (a in 1:nab) {
		for (t in 1:ntime) {
			for (d in 1:ncond){
				for (k in 1:4){
					x[a,t,d,k]=as.character(antib[a,k])				
				}				
				x[a,t,d,5]=data[(a-1)*nrep*ntime*ncond+(nrep-1)*ntime*ncond+(t-1)*ncond+d,2]	
				x[a,t,d,6]=(data[(a-1)*nrep*ntime*ncond+(1-1)*ntime*ncond+(t-1)*ncond+d,4])/(normal[a,4])
				x[a,t,d,7]=(data[(a-1)*nrep*ntime*ncond+(2-1)*ntime*ncond+(t-1)*ncond+d,4])/(normal[a,4])
				
				
			}
		}	
}

for (a in 1:nab) {
		for (t in 1:ntime) {
			for (d in 1:ncond){

					
				data_list[(a-1)*ntime*ncond+ntime*ncond+(t-1)*ncond+d,]=x[a,t,d,]			
				
			}
		}	
}
y<-as.table(x)
write.table(datalist,file="deneme.txt",quote=F,row.names=F,col.names=F)
print(antib[,1])

rm(list=ls())
