.packageName <- "zeptosensPkg"

#' @importFrom rJava .jpackage .jcall
.onLoad <- function(lib, pkg) {
  rJava::.jpackage(pkg, jars = c("causalpath.jar"))

  # Taken from xlsxjars packages What's your java version?  Need >= 1.5.0.
  jversion <- rJava::.jcall("java.lang.System", "S", "getProperty", "java.version")
  if (jversion < "1.8.0") {
    stop(paste("Your java version is ", jversion, ". Need 1.8.0 or higher.", sep = ""))
  }
}
