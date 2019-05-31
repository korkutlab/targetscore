#' Plot vocano plot of Target Score Result. As Significant Proteins (p value< Defalut or Manually set value) will show in red. The inverse log 10 of Target Score q value and Target Score calculated were shown.
#' 
#' @param TS input Target Score calculated for each protein data frame. Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param Qvalue input Target Score q value calculated for each protein data frame.Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param filename Manually set filename of Vocano Plot.
#' @param path Plot Store path. Default at working environment.
#' @return Vocano Plots with indicated filename and path.
#' @examples optimizeParameter(data=GeneExpresssion,prior=Priorindormation)
#' @concept zeptosensPkg
#' @export

getVocanoPlot=function(TS,Qvalue,filename,path=NULL){
    TS<- as.matrix(TS) 
    Padj<- as.matrix(Qvalue)
    
    if(nrow(Padj)!=nrow(TS)){
        stop("ERROR:Tag of TS and Qvalue does not match.")
    }
    
    tmpDat <- data.frame(cbind(TS,-1*log10(Padj)))
    colnames(tmpDat) <- c("TS","neglogQ")
    
    color <- ifelse(Padj>0.4,"not significant","significant")
    rownames(color) <- rownames(TS)
    tmpDat$labelnames <-  row.names(tmpDat)
    sig01 <- subset(tmpDat, tmpDat$neglogQ > -1*log10(0.4))
    siglabel <- sig01$labelnames
    tmpDat$color <- color
    
    (p <- ggplot() +
            geom_point(data=tmpDat, aes(x=TS, y=neglogQ, color=color), alpha=0.4, size=2) +
            theme_bw() +
            xlab("<TS>") + ylab("-log10 (Q-Value)") + ggtitle("")+
            scale_color_manual(name="", values=c("black", "red"))+
            geom_label_repel(data=sig01, aes(x=sig01$TS, y=sig01$neglogQ,label=siglabel), size=5)
    )
    
    plotname=paste0(path,filename,".pdf")
    ggsave(plotname,p)
    tmpDatF <- cbind(tmpDat$TS,tmpDat$neglogQ)
    colnames(tmpDatF)<-c("TS","neglogQ")
    csvname=paste0(path,filename,".csv")
    write.csv(tmpDatF,file=csvname)
}