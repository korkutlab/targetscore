#' Plot volcano plot of Target Score Result. 
#' As Significant Proteins (p value< Defalut or Manually set value) will show in red. 
#' The inverse log 10 of Target Score q value and Target Score calculated were shown.
#'
#' @param ts input Target Score calculated for each protein data frame. 
#' Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param qValue input Target Score q value calculated for each protein data frame.
#' Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param filename Manually set filename of volcano Plot.
#' @param path Plot Store path. Default at working environment.
#' @return volcano Plots with indicated filename and path.
#' @examples
#' optimizeParameter(data = GeneExpresssion, prior = Priorindormation)
#' @concept zeptosensPkg
#' @export

getVolcanoPlot <- function(ts, qValue, filename, path = NULL, sigValue = 0.4) {
  ts <- as.matrix(ts)
  pAdj <- as.matrix(qValue)

  if (nrow(pAdj) != nrow(ts)) {
    stop("ERROR:Tag of ts and qValue does not match.")
  }

  tmpDat <- data.frame(cbind(ts, -1 * log10(pAdj)))
  colnames(tmpDat) <- c("ts", "neglogQ")

  color <- ifelse(pAdj > sigValue, "not significant", "significant")
  rownames(color) <- rownames(ts)
  tmpDat$labelnames <- row.names(tmpDat)
  sig01 <- subset(tmpDat, tmpDat$neglogQ > -1 * log10(sigValue))
  siglabel <- sig01$labelnames
  tmpDat$color <- color

  (p <- ggplot() +
    geom_point(data = tmpDat, aes(x = ts, y = neglogQ, color = color), alpha = 0.4, size = 2) +
    theme_bw() +
    xlab("<ts>") + ylab("-log10 (Q-Value)") + ggtitle("") +
    scale_color_manual(name = "", values = c("black", "red")) +
    geom_label_repel(data = sig01, aes(x = sig01$ts, y = sig01$neglogQ, label = siglabel), size = 5)
  )

  plotname <- paste0(path, filename, ".pdf")
  ggsave(plotname, p)
  tmpDatF <- cbind(tmpDat$ts, tmpDat$neglogQ)
  colnames(tmpDatF) <- c("ts", "neglogQ")
  csvname <- paste0(path, filename, ".csv")
  write.csv(tmpDatF, file = csvname)
}
