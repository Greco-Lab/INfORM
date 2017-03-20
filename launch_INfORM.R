library(shiny)

#' Launch INfORM application.
#'
#' Launches the Shiny web interface for INfORM
#'
#' @import shinydashboard
#' @import DT
#' @import R.utils
#' @import png
#' @importFrom shinyjs enable disable toggle info onclick
#' @importFrom shiny runApp
#' @export
launch_INfORM <- function(){
	appDir <- system.file("INfORM-app", package = "inform")
	print(system.file("INfORM-app", package = "inform"))
	if (appDir == "") {
	stop("Unable to launch INfORM Shiny interface, please check installation!.", call. = FALSE)
	}

	shiny::runApp(appDir, display.mode = "normal")
}
