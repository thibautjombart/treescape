#'
#' Phylogenetic tree exploration
#'
#' This functions are under development. Please do not use them without contacting the author first.
#'
#' @param x
#' @param explorer
#' @param ...
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}, Michelle Kendall \email{michelle.louise.kendall@@gmail.com}, Caroline Colijn \email{c.colijn@@imperial.ac.uk}
#'
#' @export
#'
#' @import ape ade4
#' @importFrom adephylo distTips
#'
exploratree <- function(x, explorer=distTips, nf=NULL, ...){
    ## CHECKS ##
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")

    ## GET DISTANCES BETWEEN TREES ##
    ## get data.frame of all summary vectors ##
    df <- data.frame(lapply(x, function(e) as.vector(explorer(e, ...))))

    ## get pairwise Euclidean distances ##
    D <- dist(df)

    ## perform PCoA/MDS ##
    pco <- dudi.pco(D, scannf=is.null(nf), nf=nf)


    ## BUILD RESULT AND RETURN ##
    out <- list(D=D, pco=pco)
    return(out)
} # end exploratree

