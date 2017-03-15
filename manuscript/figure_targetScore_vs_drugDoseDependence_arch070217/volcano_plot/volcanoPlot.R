#volcano plot for q vs Target score
library(ggplot2)
library(ggrepel)
source("https://gist.githubusercontent.com/cannin/db5174a74349e601fbcd662f1fa2951f/raw/saveGgplotPlot.R")

#tmpDat has 2 columns FC, negLogPadj
qVal <- read.table("ovcar3_q.txt",header=T, row.names=1)
tsVal<-read.table("ovcar3_TS.txt",header=F, row.names=1)
protnames <- read.table("prot_names.txt",header=T)
negLogPadj <- -log10(qval$FDR_adjusted_p)
FC <- tsVal$V2 
    
tmpDat <- data.frame(cbind(FC,negLogPadj))
colnames(tmpDat) <- c("FC","neglogP")
color <- ifelse(qVal>0.005,"not significant","significant") 
tmpDat$labelnames <-  protnames
sig01 <- subset(tmpDat, tmpDat$neglogP > 3)
siglabel <- sig01$labelnames
#rep("significant", nrow(tmpDat))
tmpDat$color <- color
p <- ggplot() +
    geom_point(data=tmpDat, aes(x=FC, y=neglogP, color=color), alpha=0.4, size=2) +
    theme_bw() +
    xlab("TS") + ylab("-log10 (Q-Value)") + ggtitle("Ovcar3")+
    scale_color_manual(name="", values=c("black", "red"))#+
#    geom_label_repel(data=sig01, aes(x=sig01$FC, y=sig01$neglogP,label=siglabel), size=3)

ggsave("ovcar3_valcano.pdf",p)

p <- ggplot() +
  geom_point(data=tmpDat, aes(x=FC, y=negLogPadj, color=color, text=toolTip), alpha=0.4, size=2) +
  theme_bw() +
  xlab("log2 (Target score)") + ylab("-log10 Q-Value") +
  scale_color_manual(name="", values=c("red", "black"))

saveGgplotPlot(p)
