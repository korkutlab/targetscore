.packageName <- "zeptosensPkg"

#' @import rJava
.onLoad <- function(lib, pkg) {
    .jpackage(pkg, jars = c("causality.jar"))
    
    # Taken from xlsxjars packages What's your java version?  Need >= 1.5.0.
    jversion <- .jcall("java.lang.System", "S", "getProperty", "java.version")
    if (jversion < "1.8.0") {
        stop(paste("Your java version is ", jversion, ". Need 1.8.0 or higher.", sep = ""))
    }
}

#' Skip a test if on Bioconductor
#' 
#' Extension on testthat code
#' 
#' @concept paxtoolsr
#' @export
skip_on_bioc <- function() {
    if (identical(Sys.getenv("NOT_BIOC"), "true")) 
        return()
    
    message <- "On Bioconductor"
    cond <- structure(list(message = message), class = c("skip", "condition"))
    stop(cond)
} 
