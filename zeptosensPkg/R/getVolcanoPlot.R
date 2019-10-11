#' Plot volcano plot of Target Score Result.
#' As Significant Proteins (p value< Defalut or Manually set value) will show in red.
#' The inverse log 10 of Target Score q value and Target Score calculated were shown.
#'
#' @param ts input Target Score calculated for each protein data frame.
#' Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param q_value input Target Score q value calculated for each protein data frame.
#' Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param filename Manually set filename of volcano Plot.
#' @param path Plot Store path. Default at working environment.
#' @param sig_value TODO
#' 
#' @return volcano Plots with indicated filename and path.
#' 
#' @examples
#' optimizeParameter(data = GeneExpresssion, prior = Priorindormation)
#' 
#' @importFrom ggplot2 ggsave ggplot aes xlab theme_bw ggtitle xlab ylab geom_point scale_color_manual
#' @importFrom ggrepel geom_label_repel
#' @importFrom utils write.csv
#' 
#' @concept zeptosensPkg
#' @export
getVolcanoPlot <- function(ts, q_value, filename, path = NULL, sig_value = 0.4) {
  ts <- as.matrix(ts)
  p_adj <- as.matrix(q_value)

  if (nrow(p_adj) != nrow(ts)) {
    stop("ERROR:Tag of ts and q_value does not match.")
  }

  tmp_dat <- data.frame(cbind(ts, -1 * log10(p_adj)))
  colnames(tmp_dat) <- c("ts", "neglogQ")

  color <- ifelse(p_adj > sig_value, "not significant", "significant")
  rownames(color) <- rownames(ts)
  tmp_dat$labelnames <- row.names(tmp_dat)
  sig01 <- subset(tmp_dat, tmp_dat$neglogQ > -1 * log10(sig_value))
  siglabel <- sig01$labelnames
  tmp_dat$color <- color

  (p <- ggplot() +
    geom_point(data = tmp_dat, aes(x = ts, y = neglogQ, color = color), alpha = 0.4, size = 2) +
    theme_bw() +
    xlab("<ts>") + ylab("-log10 (Q-Value)") + ggtitle("") +
    scale_color_manual(name = "", values = c("black", "red")) +
    geom_label_repel(data = sig01, aes(x = sig01$ts, y = sig01$neglogQ, label = siglabel), size = 5)
  )

  plotname <- paste0(path, filename, ".pdf")
  ggplot2::ggsave(plotname, p)
  tmp_dat_f <- cbind(tmp_dat$ts, tmp_dat$neglogQ)
  colnames(tmp_dat_f) <- c("ts", "neglogQ")
  csvname <- paste0(path, filename, ".csv")
  write.csv(tmp_dat_f, file = csvname)
}
