#' Launch Famnesia
#'
#' This function launches the Famnesia Shiny application.
#'
#' @param launch.browser Whether to launch the app in a browser. Default: TRUE.
#' @param ... Further arguments passed to [shiny::runApp()], such as
#'   `quiet`, `port` and `host`
#' @return The value returned by [shiny::runApp()]
#'
#' @examples
#' if(interactive()) launchApp()
#'
#' @export
launchApp = function(launch.browser = TRUE, ...) {
  appDir = system.file("shiny", package = "famnesia", mustWork = TRUE)
  suppressPackageStartupMessages(
    shiny::runApp(appDir, launch.browser = launch.browser, ...))
}
