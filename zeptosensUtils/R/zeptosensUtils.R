.packageName <- "zeptosensUtils"

#' @importFrom rJava .jpackage
.onLoad <- function(lib, pkg) {
  # rJava::.jpackage(pkg, jars=c('causality.jar'))

  # Taken from xlsxjars packages What's your java version?  Need >= 1.5.0.
  jversion <- rJava::.jcall("java.lang.System", "S", "getProperty", "java.version")
  if (jversion < "1.5.0") {
    stop(paste("Your java version is ", jversion, ". Need 1.5.0 or higher.",
      sep = ""
    ))
  }
}
