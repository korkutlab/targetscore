#' Save a ggplot to PDF 
#' 
#' @param plt a ggplot object 
#' @param filname the plot filename
#' 
#' @export
saveGgplotPlot <- function(plt, filename = "plot.pdf") {
    pdf(filename)
    print(plt)
    dev.off()
}