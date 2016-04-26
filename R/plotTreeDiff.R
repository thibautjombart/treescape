#' Plot tree differences
#'
#' Highlight the topologicial differences between two trees, plotted side by side. 
#' This function is useful for comparing representative "median" trees - see \code{\link{medTree}}.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom combinat combn2
#'
#' @param tr1 an object of the class \code{phylo}: the first tree to plot
#' @param tr2 an object of the class \code{phylo}: the second tree to plot
#' @param vec1 an optional input, the result of \code{treeVec(tr1, lambda=0)}, to speed up the computation
#' @param vec2 an optional input, the result of \code{treeVec(tr2, lambda=0)}, to speed up the computation
#' @param baseCol the colour used for tips with identical ancestry in the two trees. Defaults to "grey".
#' @param col1 the first colour used to define the colour spectrum for tips with differences. This colour will be used for tips with minor differences. Defaults to "peachpuff".
#' @param col2 the second colour used to define the colour spectrum for tips with differences. This colour will be used for tips with major differences. Defaults to "red2".
#' @param ... further arguments passed to \code{\link{plot.phylo}}
#'
#' @return
#' A plot of the two trees side by side. Tips are coloured in the following way:
#' \itemize{
#' \item if each ancestor of a tip in tree 1 occurs in tree 2 with the same partition of tip descendants, then the tip is coloured grey (or supplied "baseCol")
#' \item if not, the tip gets coloured pale orange to red on a scale according to how many differences there are amongst its most recent common ancestors with other tips. The colour spectrum can be changed according to preference.
#' }
#' 
#' @seealso \code{\link{medTree}}
#'
#' @examples
#' ## simple example on trees with five tips:
#' tr1 <- read.tree(text="((A:1,B:1):1,((C:1,D:1):1,E:1):1):1;")
#' tr2 <- read.tree(text="((A:1,B:1):1,(C:1,(D:1,E:1):1):1):1;")
#' plotTreeDiff(tr1,tr2)
#' 
#' ## example on larger woodmice trees
#' data(woodmiceTrees)
#' plotTreeDiff(woodmiceTrees[[1]],woodmiceTrees[[2]])
#' ## change aesthetics:
#' plotTreeDiff(woodmiceTrees[[1]],woodmiceTrees[[2]], 
#'    baseCol="grey2", col1="cyan", col2="navy", 
#'    edge.width=2, type="radial", cex=0.5, font=2)
#' 
#' @export 
plotTreeDiff <- function(tr1,tr2,vec1=NULL,vec2=NULL,baseCol="grey",col1="peachpuff",col2="red2",...) {
  
  l <- length(tr1$tip.label)
  lchoose2 <- l*(l-1)/2
  if( l != length(tr2$tip.label)) stop("Trees must have the same number of tips")
  
  if(setequal(tr1$tip.label,tr2$tip.label) == FALSE) stop("Trees must have the same tip label sets")
  
  # get vec1
  if (is.null(vec1)) {
    vec1 <- treeVec(tr1) # emphasise.tips, emphasise.weight)
  }
  if (is.null(vec2)) {
    vec2 <- treeVec(tr2) # emphasise.tips, emphasise.weight)
  }
  
  # trim pendant edge entries from vectors
  vec1 <- vec1[1:lchoose2] # emphasise.tips, emphasise.weight)
  vec2 <- vec2[1:lchoose2]
  
  # find the positions where the vectors are different
  vecDiff <- (vec1!=vec2)
  
  # combine the tip labels (in order) and whether the vectors are different
  treedf <- as.data.frame(cbind(combn2(tr1$tip.label[order(tr1$tip.label)]),vecDiff))
  
  # list the tips which appear in "TRUE" rows
  tipDiffs <- c(as.character(treedf[,1][vecDiff]),as.character(treedf[,2][vecDiff]))
  
  # find the number of times each tip appears, in the order that the tip labels are read
  # add one because "zero" won't be defined in the colour palette
  tipSignificance1 <- sapply(tr1$tip.label, function(x)
    length(which(tipDiffs==x))+1)
  tipSignificance2 <- sapply(tr2$tip.label, function(x)
    length(which(tipDiffs==x))+1)
  
  numCols <- max(tipSignificance1) -1 
  colfunc<-colorRampPalette(c(col1,col2))
  pal <- c(baseCol,colfunc(numCols))
  
  tipCols1 <- pal[tipSignificance1]
  tipCols2 <- pal[tipSignificance2]
  
  # plot
  layout(matrix(1:2, 1, 2))
  plot.phylo(tr1, tip.color=tipCols1, no.margin=TRUE, ...)
  plot.phylo(tr2, tip.color=tipCols2, no.margin=TRUE, ...)
  
}