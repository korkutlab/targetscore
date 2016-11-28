#' Save a ggplot to PNG/PDF 
#' 
#' @param plt a ggplot object 
#' @param type either "png" or "pdf"
#' @param filename the plot filename
#' @param ... additional parameters for either png() or pdf() 
#' 
#' @export
saveGgplotPlot <- function(plt, type=c("png", "pdf"), filename="plot", ...) {
    filename <- paste0(filename, type)
    
    if(type == "png") {
        png(filename, ...)
    } else if (type == "pdf") {   
        pdf(filename, ...)
    }
    
    print(plt)
    dev.off()
}