#'
#' Phylogenetic tree exploration
#'
#' Compares phylogenetic trees and maps them into a small number of dimensions for easy visualisation and identification of clusters.
#'
#' @param x an object of the class multiPhylo
#' @param method a function outputting the summary of a tree (phylo object) in the form of a vector
#' @param nf the number of principal components to retain
#' @param ... further arguments to be passed to \code{method}
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}, Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @export
#'
#' @import ape
#' @importFrom ade4 dudi.pco
#' @importFrom adephylo distTips
#' @importFrom stats dist
#'
#' @examples
#'
#' ## generate list of trees
#' x <- rmtree(10, 20)
#' names(x) <- paste("tree", 1:10, sep = "")
#'
#' ## use treescape
#' res <- treescape(x, nf=3)
#' table.paint(as.matrix(res$D))
#' scatter(res$pco)
#'
#' data(woodmiceTrees)
#' woodmiceDists <- treescape(woodmiceTrees,nf=3)
#' plot(woodmiceDists$pco$li[,1],woodmiceDists$pco$li[,2])
#' woodmicedf <- woodmiceDists$pco$li
#' if(require(ggplot2)){
#' woodmiceplot <- ggplot(woodmicedf, aes(x=A1, y=A2)) # create plot
#' woodmiceplot + geom_density2d(colour="gray80") + # contour lines
#' geom_point(size=6, shape=1, colour="gray50") + # grey edges
#' geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
#' xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
#' }
#' \dontrun{
#' if(require(rgl)){
#' plot3d(woodmicedf[,1], woodmicedf[,2], woodmicedf[,3], type="s", size=1.5,
#' col="navy", alpha=0.5, xlab="", ylab="", zlab="")
#' }
#' }
treescape <- function(x, method=treeVec, nf=NULL, ...){
    ## CHECKS ##
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")
    if(is.null(names(x))) names(x) <- 1:length(x)
    lab <- names(x)


    ## GET DISTANCES BETWEEN TREES ##
    ## get data.frame of all summary vectors ##
    df <- t(data.frame(lapply(x, function(e) as.vector(method(e, ...)))))

    ## get pairwise Euclidean distances ##
    D <- dist(df)
    attr(D,"Labels") <- lab # restore labels

    ## perform PCoA/MDS ##
    pco <- dudi.pco(D, scannf=is.null(nf), nf=nf)


    ## BUILD RESULT AND RETURN ##
    out <- list(D=D, pco=pco)
    return(out)
} # end treescape
