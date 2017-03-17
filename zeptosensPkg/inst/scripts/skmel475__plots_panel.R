#plot proteomics data
library(Hmisc)
library(graphics)
data <- read.table("skmel_475_a_n_tp.txt")
nab=72
ntime=3
ncond=4
pdf("panel1.pdf")

par(mfrow=c(4,4),mar = rep(2, 4))
for(i in 1:16) {
	time = data[1:ntime,5]
readout = data[((i-1)*ntime*ncond+(ncond-3)*ntime+1):((i-1)*ntime*ncond+(ncond-3)*ntime+ntime),8]
err=data[((i-1)*ntime*ncond+(ncond-3)*ntime+1):((i-1)*ntime*ncond+(ncond-3)*ntime+ntime),9]
readout2 = data[((i-1)*ntime*ncond+(ncond-2)*ntime+1):((i-1)*ntime*ncond+(ncond-2)*ntime+ntime),8]
err2=data[((i-1)*ntime*ncond+(ncond-3)*ntime+1):((i-1)*ntime*ncond+(ncond-3)*ntime+ntime),9]
readout3 = data[((i-1)*ntime*ncond+(ncond-1)*ntime+1):((i-1)*ntime*ncond+(ncond-1)*ntime+ntime),8]
err3=data[((i-1)*ntime*ncond+(ncond-3)*ntime+1):((i-1)*ntime*ncond+(ncond-3)*ntime+ntime),9]
readout
opts = c("l")
title=data[(i-1)*ntime*ncond+(ncond-3)*ntime+1,7]
#extension=".pdf"
#pdf(paste(toString(title, width=NULL),extension))
#plot(time,readout,type=opts[1])
errbar(time,readout,type=opts[1],readout+err,readout-err,xlim=c(0,24),ylim=c(0,3),col=rgb(1,0,0))
par(new=T)
errbar(time,readout2,type=opts[1],readout2+err2,readout2-err2,xlim=c(0,24),ylim=c(0,3),col=rgb(0,1,0),ann=FALSE)
par(new=T)
errbar(time,readout3,type=opts[1],readout3+err3,readout3-err3,xlim=c(0,24),ylim=c(0,3),col=rgb(0,0,1),xlab="time")

title(main=title,cex=0.3)

#dev.off()
}
dev.off()
