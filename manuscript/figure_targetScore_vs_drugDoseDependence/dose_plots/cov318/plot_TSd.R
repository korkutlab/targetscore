library(ggplot2)
#input TS and stdev 
source("https://gist.githubusercontent.com/cannin/db5174a74349e601fbcd662f1fa2951f/raw/saveGgplotPlot.R")
data1 <- read.table(file="cov318_TS_d_p1.txt",header=T)
data0 <- read.table(file="cov318_TS_d_p0.txt",header=T) 
err1 <-  read.table(file="cov318_TS_stdev_p1.txt",header=T)
err0 <-  read.table(file="cov318_TS_stdev_p0.txt",header=T)
prot_list <- read.table(file="prot_list",header=T)

#extract nodes of interest from match btwn prot_list and data
plot_data1 <- data1[,colnames(data1) %in% prot_list$proteins]
plot_data0 <- data0[,colnames(data0) %in% prot_list$proteins]
plot_err1 <- err1[,colnames(err1) %in% prot_list$proteins]
plot_err0 <- err0[,colnames(err0) %in% prot_list$proteins]
#revert the order of rows to have ascending doses
plot_data1 <- plot_data1[nrow(plot_data1):1,]
plot_data0 <- plot_data0[nrow(plot_data0):1,]
plot_err1 <- plot_err1[nrow(plot_err1):1,]
plot_err0 <- plot_err0[nrow(plot_err0):1,]
#plot
dose <- c(0.01,0.02,0.04,0.08,0.16,0.31,0.63,1.25,2.5,5)
#plotdata
plot_data1x <- cbind(dose,plot_data1)
plot_data0x <- cbind(dose,plot_data0)
plot_data1x$sup <- rep("p1",each=10)
plot_data0x$sup <- rep("p0",each=10)
plot_data <- rbind(plot_data0x,plot_data1x)

#ploterror
plot_err1x <- cbind(dose,plot_err1)
plot_err0x <- cbind(dose,plot_err0)
plot_err1x$sup <- rep("p1",each=10)
plot_err0x$sup <- rep("p0",each=10)
plot_err <- rbind(plot_err0x,plot_err1x)
for (i in 1:6){
#  pdf("a.pdf")
plt <-  ggplot(data=plot_data, aes(x=dose, y=plot_data[,i+1], group=sup)) +
    geom_line()+
    geom_point(aes(shape=sup))+
    geom_line(aes(color=sup))+
    geom_point(aes(color=sup))+
    scale_x_log10()+
    theme_minimal()+
    ylab("TSd")+
    xlab("Dose (uM)")+
    geom_errorbar(aes(ymin=plot_data[,i+1]-0.57*plot_err[,i+1], ymax=plot_data[,i+1]+0.57*plot_err[,i+1]), width=.1)
ggsave(paste0(colnames(plot_data)[i+1],".pdf"),plot=plt)
  }



