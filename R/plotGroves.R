#'
#' Scatterplot of groups of trees
#'
#' This function displays the scatterplot of the Multidimensional
#' Scaling (MDS) output by treescape, superimposing group information
#' (derived by \code{\link{findGroves}}) using colors.
#'
#' This function relies on \code{\link[adegraphics]{s.class}}
#' from the adegraphics package.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @importFrom adegraphics s.class
#' @importFrom adegraphics s.label
#' @importFrom adegraphics s1d.barchart
#' @importFrom adegraphics insert
#' @importFrom adegenet funky
#' @importFrom adegenet bluepal
#' @importFrom adegenet transp
#'
#' @param x a list returned by \code{\link{findGroves}} or a MDS with class \code{dudi}
#' @param groups a factor defining groups of trees
#' @param xax a number indicating which principal component to be used as 'x' axis
#' @param yax a number indicating which principal component to be used as 'y' axis
#' @param type a character string indicating which type of graph to use
#' @param col.pal a color palette to be used for the groups
#' @param bg the background color
#' @param lab.show a logical indicating whether labels should be displayed
#' @param lab.col a color for the labels
#' @param lab.cex the size of the labels
#' @param lab.optim a logical indicating whether label positions should be optimized to avoid overlap; better display but time-consuming for large datasets
#' @param point.cex the size of the points
#' @param scree.pal a color palette for the screeplot
#' @param scree.size a size factor for the screeplot, between 0 and 1
#' @param scree.posi either a character string or xy coordinates indicating the position of the screeplot.
#' @param ... further arguments passed to \code{\link{s.class}}
#'
#' @return
#' An adegraphics object (class: ADEgS)
#'
#' @seealso
#' \code{\link[adegraphics]{s.class}}
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' if(require("adegenet") && require("adegraphics")){
#' ## load data
#' data(woodmiceTrees)
#'
#' ## run findGroves: treescape+clustering
#' res <- findGroves(woodmiceTrees, nf=5, nclust=6)
#'
#' ## basic plot
#' plotGroves(res)
#'
#' ## adding labels
#' plotGroves(res, lab.show=TRUE)
#'
#' ## customizing
#' plotGroves(res, lab.show=TRUE,
#' bg="black", lab.col="white", scree.size=.35)
#'
#' ## customizing
#' plotGroves(res, type="ellipse", lab.show=TRUE,
#' lab.optim=FALSE, scree.size=.35)
#'
#' ## example with no group information
#' plotGroves(res$treescape$pco)
#'
#' ## adding labels
#' plotGroves(res$treescape$pco, lab.show=TRUE, lab.cex=2)
#'
#' }
#' }
#'
plotGroves <- function(x, groups=NULL, xax=1, yax=2,
                        type=c("chull","ellipse"), col.pal=funky, bg="white",
                        lab.show=FALSE, lab.col="black", lab.cex=1, lab.optim=TRUE,
                        point.cex=1, scree.pal=NULL, scree.size=.2,
                        scree.posi=c(.02,.02), ...){
    ## HANDLE ARGUMENTS ##
    ## checks
    type <- match.arg(type)
    if(is.null(scree.pal)) scree.pal <- function(n) rev(bluepal(n))

    ## x is a list returned by findGroves
    if(is.list(x) && !is.data.frame(x) && !inherits(x,"dudi")){
        if(is.null(x$treescape)) stop("if x is a list, it should contain a slot $treescape")
        groups <- x$groups
        x <- x$treescape$pco
    }

    ## x is a dudi object
    if(inherits(x,"dudi")){
        eig <- x$eig
        x <- x$li
    }

    ## groups missing - just s.label
    if(is.null(groups)) {
        ## with labels
        if(lab.show){
            out <- s.label(x, xax=xax, yax=yax,
                           plabels=list(optim=lab.optim, col=lab.col, cex=lab.cex),
                           pppoints=list(cex=point.cex),
                           pbackground.col=bg,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        } else {
            ## just points
            out <- s.label(x, xax=xax, yax=yax,
                           plabels=list(optim=FALSE,cex=0),
                           ppoints=list(cex=point.cex, col=lab.col),
                           pbackground.col=bg,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        }
    } else {
        ## if groups are provided
        if(!is.factor(groups)) groups <- factor(groups)
        n.lev <- length(levels(groups))


        ## MAKE GRAPH ##
        ## base scatterplot
        if(type=="chull"){
            out <- s.class(x, xax=xax, yax=yax, fac=groups, col=col.pal(n.lev),
                           ellipseSize=0, chullSize=1,
                           pbackground.col=bg,
                           ppoints.cex=point.cex,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        }
        if(type=="ellipse"){
            out <- s.class(x, xax=xax, yax=yax, fac=groups, col=col.pal(n.lev),
                           ellipseSize=1,
                           pbackground.col=bg,
                           ppoints.cex=point.cex,
                           pgrid.text.col=lab.col, plot=FALSE, ...)
        }

        ## add labels
        if(lab.show){
            out <- out + s.label(x, plabel.optim=lab.optim, plabel.col=lab.col,
                                 ppoints.cex=0, plabels.cex=lab.cex)
        }
    }
    ## add inset
    if(!is.null(scree.posi[1]) && !is.na(scree.posi[1]) && scree.posi[1]!="none"){
        screeplot <- s1d.barchart(c(rep(0,3),eig), p1d.horizontal=FALSE, ppolygons.col=scree.pal(length(eig)),
                                  pbackground=list(col=transp("white"), box=TRUE),
                                  layout.width=list(left.padding=2),
                                  pgrid.draw=FALSE, plot=FALSE)
        out <- insert(screeplot, out, posi=scree.posi, ratio=scree.size, plot=FALSE)

    }


    ## RETURN ##
    return(out)
} # end plotGroves
