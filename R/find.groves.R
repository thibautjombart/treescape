#'
#' Identify clusters of similar trees
#'
#' This function uses hierarchical clustering on principal components output by \code{\link{treescape}} to identify groups of similar trees. Clustering relies on \code{\link{hclust}}, using Ward's method by default.
#'
#' @param x an object of the class multiPhylo
#' @param method a function outputting the summary of a tree (phylo object) in the form of a vector
#' @param nf the number of principal components to retain
#' @param clustering a character string indicating the clustering method to be used; defaults to Ward's method; see argument \code{method} in \code{?hclust} for more details.
#' @param nclust an integer indicating the number of clusters to find; if not provided, an interactive process based on cutoff threshold selection is used.
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
#' ## run find.groves: treescape+clustering
#' res <- find.groves(woodmiceTrees, nf=5, nclust=6)
#'
#' ## plot results on first 2 axes
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





#' Scatterplot of trees
#'
#' This function displays the scatterplot of the Multidimensional Scaling (MDS) output by treescape, superimposing group information (derived by \code{\link{find.groves}} using colors.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @importFrom adegraphics s.class
#' @importFrom adegraphics s1d.barchart
#' @importFrom adegraphics insert
#' @importFrom adegenet funky
#' @importFrom adegenet bluepal
#' @importFrom adegenet transp
#'
plot.groves <- function(x, groups=NULL, xax=1, yax=2,
                        type=c("chull","ellipse"), col.pal=funky,
                        show.labels=FALSE,
                        bg=transp("black"), lab.col="white",
                        scree.pal=NULL, scree.size=.2,
                        scree.posi=c(.02,.02), ...){
    ## HANDLE ARGUMENTS ##
    ## checks
    type <- match.arg(type)
    if(is.null(scree.pal)) scree.pal <- function(n) rev(bluepal(n))

    ## x is a list returned by find.groves
    if(is.list(x) && !is.data.frame(x) && !inherits(x,"dudi")){
        if(is.null(x$groups)) stop("if x is a list, it should contain a slot $groups")
        if(is.null(x$treescape)) stop("if x is a list, it should contain a slot $treescape")
        groups <- x$groups
        x <- x$treescape$pco
    }

    ## x is a dudi object
    if(inherits(x,"dudi")){
        eig <- x$eig
        x <- x$li
    }

    ## groups
    if(is.null(groups)) stop("group information missing; try running find.groves first")
    if(!is.factor(groups)) groups <- factor(groups)
    n.lev <- length(levels(groups))


    ## MAKE GRAPH ##
    ## base scatterplot
    if(type=="chull"){
        out <- s.class(x, xax=xax, yax=yax, fac=groups, col=col.pal(n.lev),
                       ellipse=0, chullSize=1, pbackground.col=bg, plot=FALSE)
    }
    if(type=="ellipse"){
        out <- s.class(x, xax=xax, yax=yax, fac=groups, col=col.pal(n.lev),
                       pbackground.col=bg, ellipse=1, plot=FALSE)
    }

    ## add labels
    if(show.labels){
        out <- out + s.label(x$li, plabel.optim=TRUE, plabel.col=lab.col, ppoints.cex=0)
    }

    ## add inset
    if(!is.null(scree.posi[1]) && !is.na(scree.posi[1])){
        screeplot <- s1d.barchart(eig, p1d.horizontal=FALSE, ppolygons.col=scree.pal(length(eig)),
                                  pbackground.col=transp("white"), pgrid.draw=FALSE, plot=FALSE)
        out <- insert(screeplot, out, posi=scree.posi, ratio=scree.size, plot=FALSE)

    }


    ## RETURN ##
    return(out)
} # end plot.groves
