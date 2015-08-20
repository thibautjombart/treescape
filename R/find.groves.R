#'
#' Identify clusters of similar trees
#'
#' This function uses hierarchical clustering on principal components output by \code{\link{treescape}} to identify groups of similar trees. Clustering relies on \code{\link{hclust}}, using Ward's method by default.
#'
#' @param x an object of the class multiPhylo
#' @param method a function outputting the summary of a tree (phylo object) in the form of a vector
#' @param nf the number of principal components to retain
#' @param clustering a character string indicating the clustering method to be used; see argument \code{method} in \code{?hclust} for more details.
#' @param ... further arguments to be passed to \code{treescape}
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @import ape
#' @importFrom stats hclust
#'
#' @return
#' A list containing:
#' \itemize{
#'  \item groups: a factor defining groups of trees
#'  \item treescape: the output of treescape
#' }
#'
#' @examples
#'
#' if(require("adegenet")){
#' ## load data
#' data(woodmiceTrees)
#'
#' ## use treescape
#' res <- find.groves(woodmiceTrees, nf=5, nclust=6)
#' PCs <- res$treescape$pco$li
#' s.class(PCs, fac=res$groups, col=funky(6))
#' }
#'
find.groves <- function(x, method=tree.vec, nf=NULL, clustering="ward.D2",
                        nclust=NULL, ...){
    ## CHECKS ##
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")

    ## GET OUTPUT OF TREESCAPE ##
    res <- treescape(x, method=method, nf=nf, ...)

    ## GET CLUSTERS ##
    ## hierharchical clustering
    clust <- hclust(dist(res$pco$li), method=clustering)

    ## select nclust interactively if needed
    if(is.null(nclust)){
        ans <- NA
        continue <- TRUE
        while(is.na(ans) || continue){
            plot(clust)
            cat("\nPlease define a cutoff point: ")
            ans <- as.double(readLines(n = 1))
            abline(h=ans, col="royalblue", lty=2)
            cat("\nAre you happy with this choice (y/n): ")
            continue <- as.character(readLines(n = 1))!="y"
        }
        grp <- cutree(clust, h=ans)
    } else {
        ## cut tree
        grp <- cutree(clust, k=nclust)
    }

    ## BUILD RESULT AND RETURN ##
    names(grp) <- names(x)
    out <- list(groups=factor(grp), treescape=res)

    return(out)
} # end find.groves
