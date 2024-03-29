
#' Library search GUI
#' @export
runGUI <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "MergeION")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}