#' Run Shiny App
#' 
#' @param port the port number to use (Default: 3838)
#' 
#' @examples
#' \dontrun{
#'     runShinyApp()
#' }
#' 
#' @concept kircImmuneProject
#' @export
runShinyApp <- function (port=3838) {
    shiny::runApp(system.file('shinyApp', package='kircImmuneProject'), port=port)	
}
