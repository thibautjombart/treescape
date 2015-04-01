#'
#' Phylogenetic tree exploration
#'
#' This functions are under development. Please do not use them without contacting the author first.
#'
#' @param x an object of the class \code{\link[ape]{multiPhylo}}
#' @param method a function outputting the summary of a tree (phylo object) in the form of a vector
#' @param nf the number of principal components to retain
#' @param ... further arguments to be passed to \code{method}
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}, Michelle Kendall \email{michelle.louise.kendall@@gmail.com}, Caroline Colijn \email{c.colijn@@imperial.ac.uk}
#'
#' @export
#'
#' @import ape ade4
#' @importFrom adephylo distTips
#'
#' @examples
#'
#' ## generate list of trees
#' x <- rmtree(10, 20)
#' names(x) <- paste("tree", 1:10, sep = "")
#'
#' ## use exploratree
#' res <- exploratree(x, nf=3)
#' table.paint(as.matrix(res$D))
#' scatter(res$pco)
#'
exploratree <- function(x, method=distTips, nf=NULL, ...){
    ## CHECKS ##
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")

    ## GET DISTANCES BETWEEN TREES ##
    ## get data.frame of all summary vectors ##
    df <- t(data.frame(lapply(x, function(e) as.vector(method(e, ...)))))

    ## get pairwise Euclidean distances ##
    D <- dist(df)

    ## perform PCoA/MDS ##
    pco <- dudi.pco(D, scannf=is.null(nf), nf=nf)


    ## BUILD RESULT AND RETURN ##
    out <- list(D=D, pco=pco)
    return(out)
} # end exploratree

