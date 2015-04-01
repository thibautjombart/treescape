#'
#' Web-based tree explorer
#'
#' This function opens up an application in a web browser for an interactive exploration of the diversity in a set of trees. This function is under development. Please do not use without contacting the author.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#' @import shiny
#' @import adephylo
#' @import ade4
#' @import adegenet
#' @importFrom utils packageDescription
#'
exploratreeServer <- function(){
    runApp(system.file("shiny",package="exploratree"))
    return(invisible())
}





#'
#' Auxiliary functions
#'
#' These functions are not supposed to be used by the user.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @importFrom adegenet .readExt
#'
.render.server.info <- function(){
    renderPrint(
            {
                cat("\n== R version ==\n")
                print(R.version)

                cat("\n== Date ==\n")
                print(date())

                cat("\n== exploratree version ==\n")
                print(packageDescription("exploratree", fields=c("Package", "Version", "Date", "Built")))

                cat("\n== shiny version ==\n")
                print(packageDescription("shiny", fields=c("Package", "Version", "Date", "Built")))

                cat("\n== attached packages ==\n")
                print(search())
            }
                ) # end renderPrint
} # end .render.server.info
