#' Plot volcano plot of Target Score Result.
#' As Significant Proteins (p value< Defalut or Manually set value) will show in red.
#' The inverse log 10 of Target Score q value and Target Score calculated were shown.
#'
#' @param ts input Target Score calculated for each antibody data frame.
#' Gene in columns and samples in row. With colnames as gene tags and rownames as sample tags.
#' @param q_value input Target Score q value calculated for each antibody data frame.
#' Gene in columns and samples in row. With colnames as gene tags and rownames as sample tags.
#' @param filename Manually set filename of volcano plot.
#' @param path Plot Store path. Default at working environment.
#' @param sig_value  Manually set significant cut-off value for log10(q_value). (Default: 0.4)
#' @param sig_TS Manually set significant cut-off value for calculated Target Score Value. (Default: 0.5)
#' @param x_min Manually set minimum value for x lab. (Default: -2)
#' @param x_max Manually set maximum value for x lab. (Default: 2)
#' @param include_labels a boolean whether to point labels
#' @param save_output a boolean whether to save plots and point data to file
#' @param label_names a vector of strings to override the existing labels 
#' @param include_cutoffs a boolean whether cutoffs should be shown with additional lines
#'
#' @return volcano plots as ggplot object; plots and data maybe be saved as well
#'
#' @importFrom ggplot2 ggsave ggplot aes xlab theme_bw ggtitle xlab ylab geom_point scale_color_manual aes_string geom_hline geom_vline
#' @importFrom ggrepel geom_label_repel
#' @importFrom utils write.csv
#'
#' @concept targetscore
#' @export
get_volcano_plot <- function(ts, 
                             q_value,
                             filename,
                             path = getwd(),
                             sig_value = 0.4,
                             sig_TS = 0.5,
                             include_labels = TRUE,
                             save_output = TRUE,
                             x_min = -2,
                             x_max = 2,
                             label_names = NULL,
                             title = "",
                             include_cutoffs = FALSE
                             ) {
  ts <- as.matrix(ts)
  p_adj <- as.matrix(q_value)

  if (nrow(p_adj) != nrow(ts)) {
    stop("ERROR: Size of ts and q_value does not match.")
  }

  tmp_dat <- data.frame(ts=ts, neglogQ=-1*log10(p_adj))
  #colnames(tmp_dat) <- c("ts", "neglogQ")

  color <- ifelse(p_adj > sig_value, "Not Significant", "Significant")
  rownames(color) <- rownames(ts)
  
  if(!is.null(label_names)) {
    tmp_dat$label_name <- label_names
  } else {
    tmp_dat$label_name <- row.names(tmp_dat)
  }
  
  sig_neglog <- -1 * log10(sig_value)
  sig01 <- subset(tmp_dat, tmp_dat$neglogQ > sig_neglog)
  sig001<- subset(sig01, abs(sig01$ts) > sig_TS)
  siglabel <- sig001$label_name
  tmp_dat$color <- color
  tmp_dat <- tmp_dat[complete.cases(tmp_dat), ]
  
  # label = "label_name", 
  p <- ggplot() +
    geom_point(data = tmp_dat, aes_string(x = "ts", y = "neglogQ", color = "color"), alpha = 0.4, size = 2) +
    xlab("<TS>") +
    ylab("-log10 (Q-Value)") +
    scale_color_manual(name = "", values = c("black", "red")) +
    theme_bw()
  
  if(include_cutoffs) {
    p <- p + 
      geom_hline(yintercept = sig_neglog, color = "red") + 
      geom_vline(xintercept = sig_TS, color = "red") + 
      geom_vline(xintercept = -1*sig_TS, color = "red")
  }

  if(include_labels) {
    p <- p + geom_label_repel(data = sig001, aes_string(x = "ts", y = "neglogQ", label = "label_name"), size = 5)
  }
  
  if(title != "") {
    p <- p + ggtitle(title)
  }
  
  if(save_output) {
    plotname <- file.path(path, paste0(filename, ".pdf"))
    ggplot2::ggsave(plotname, p)
    
    tmp_dat_f <- cbind(tmp_dat$ts, tmp_dat$neglogQ)
    colnames(tmp_dat_f) <- c("ts", "neglogQ")
    csvname <- file.path(path, paste0(filename, ".csv"))
    write.csv(tmp_dat_f, file = csvname, row.names = FALSE, quote = FALSE)
  }

  return(p)
}
