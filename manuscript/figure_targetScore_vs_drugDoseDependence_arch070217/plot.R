# Basic plot
# plot(ovcar3_TS_all$X4, -ovcar3_TS_all$X3)
# cor(ovcar3_TS_all$X4, -ovcar3_TS_all$X3, use = "pairwise.complete.obs", method = "spearman")
# abline(lm(X3 ~ X4, ovcar3_TS_all))

library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(dplyr)

source("https://gist.githubusercontent.com/cannin/db5174a74349e601fbcd662f1fa2951f/raw/saveGgplotPlot.R")

ptCol <- brewer.pal(8, "Set1")[2]
ptCol <- "black"
folder <- "../manuscript/figure_targetScore_vs_drugDoseDependence/"
cellLines <- c("cov318", "ovcar4", "ovcar3", "igrov1")

abOfInterest <- c("ACC_pS79.R.V_GBL9026172", "ACC1.R.C_GBL9026173")

for(i in 1:length(cellLines)) {
    #i <- 2
    
    #sprFile <- paste0(folder, cellLines[i], "_spR_TS.txt")
    #tsFile <- paste0(folder, cellLines[i], "_TS.txt")
    datFile <- paste0(folder, cellLines[i], ".txt")
    #spr <- read.table(sprFile, sep="\t", header=TRUE, na.strings = "NaN", stringsAsFactors = FALSE)
    #ts <- read.table(tsFile, sep="\t", header=TRUE, na.strings = "NaN", stringsAsFactors = FALSE)
    dat <- read.table(datFile, sep="\t", header=TRUE, na.strings = "NaN", stringsAsFactors = FALSE)
    #dat <- merge(spr, ts)
    abOfInterest <- sample(dat$abname, 10)
    
    x1 <- strsplit(dat$abname, "\\.")
    x2 <- lapply(x1, function(x) {
        idx <- which(x %in% c("R", "M", "G"))
        #cat(x, "\n")
        paste(x[1:idx-1], collapse="_")
    })
    dat$displayName <- unlist(x2)

    dat$negSpCorr <- -dat$spCorr
    p <- ggplot(dat, aes(negSpCorr, TS))
    p <- p + geom_point(colour = ptCol) + ggtitle("Target Score") + xlab("Dose Dependent Response") + ylab("Target Score") + geom_hline(yintercept=0) + geom_vline(xintercept=0) + theme_bw()
    #p
    #p <- p + annotate("text", label="Shane Victorino", x=0.5, y=10, size=3) + annotate("segment", x=0.7, y=10, xend=0.8, yend=-12, size=0.5, arrow=arrow(length=unit(.2, "cm")))
    
    filterDat <- filter(dat, negSpCorr > 0.8 & TS >= quantile(d1$TS, 0.9))
    p <- p + geom_point(data=filterDat, color="red")
    #p <- p + geom_label_repel(data=filter(dat, dat$abname %in% abOfInterest), aes(label=displayName))
    p <- p + geom_label_repel(data=filterDat, aes(label=displayName))

    #p <- p + geom_text_repel(aes(negSpCorr, TS, label=abname))
    
    plotFilename <- paste0(folder, cellLines[i], ".pdf")
    saveGgplotPlot(p, plotFilename)
    #p
}


