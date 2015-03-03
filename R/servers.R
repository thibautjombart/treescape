#'
#' Web-based tree explorer
#'
#' This function opens up an application in a web browser for an interactive exploration of the diversity in a set of trees. This function is under development. Please do not use without contacting the author.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#' @import shiny
#' @importFrom utils packageDescription
#'
#' @param x a multiPhylo object
exploratreeServer <- function(x){
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")

    runApp(system.file("shiny",package="exploratree"))

    return(invisible())
}





#######################
## .render.server.info
#######################
## INVISIBLE FUNCTION RENDERING SERVER INFO ##
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
