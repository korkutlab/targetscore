library(shiny)

tags <- shiny::tags

# Show loading modal when the webpage is ready (different from when the app is ready)
# and prevent user from removing it.
onload_js <- '
$(document).ready(function() {
  $("#loadingModal").modal({backdrop: "static", keyboard: false, show: false});
});
'

# Simple load modal using Bootstrap
# NOTE: Shiny uses Bootstrap 3.3.1
loading_modal <- function(message) {
  modal <- tags$div(
    class = "modal",
    id = "loadingModal",
    tags$div(
      class = "modal-dialog modal-sm",
      tags$div(
        class = "modal-content",
        tags$div(
          class = "modal-header",
          tags$h4(
            class = "modal-title",
            message
          )
        ),
        tags$div(
          class = "modal-body",
          tags$div(
            class = "progress progress-striped active",
            style = "margin-bottom: 0;",
            tags$div(
              class = "progress-bar progress-bar-info",
              style = "width: 100%"
            )
          )
        )
      )
    )
  )

  tags$html(modal)
}
