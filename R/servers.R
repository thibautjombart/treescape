#'
#' Web-based tree explorer
#'
#' This function opens up an application in a web browser for an interactive exploration of the diversity in a set of trees. 
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape

#' @import shiny
#' @import rglwidget
#' @importFrom ade4 dudi.pco
#' @importFrom adephylo distTips
#' @importFrom adegraphics s.class 
#' @importFrom adegraphics s.label 
#' @importFrom adegraphics s1d.barchart
#' @importFrom phangorn Children
#' @importFrom phangorn Descendants
#' @importFrom utils packageDescription
#' @importFrom RLumShiny jscolorInput
#' @importFrom shinyBS bsTooltip
#' @importFrom shinyRGL webGLOutput
#' 
#' 
#' @export
treescapeServer <- function(){
    ## RUN APP
    runApp(system.file("shiny",package="treescape"))
    return(invisible())
}





#'
#' Auxiliary functions
#'
#' These functions are not supposed to be used by the user.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#'
#' @importFrom adegenet .readExt
#'
#' @export
.render.server.info <- function(){
    renderPrint(
            {
                cat("\n== R version ==\n")
                print(R.version)

                cat("\n== Date ==\n")
                print(date())

                cat("\n== treescape version ==\n")
                print(packageDescription("treescape", fields=c("Package", "Version", "Date", "Built")))

                cat("\n== shiny version ==\n")
                print(packageDescription("shiny", fields=c("Package", "Version", "Date", "Built")))

                cat("\n== attached packages ==\n")
                print(search())
            }
                ) # end renderPrint
} # end .render.server.info
