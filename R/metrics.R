
#' Linear MRCA function
#'
#' Function to make the MRCA matrix of a tree, where entry (i,j) gives the MRCA of tips i and j.
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @export
#'
#' @param tree an object of the class \code{phylo} 
#' @param k (optional) number of tips in tree, for faster computation
#'
#' @importFrom phangorn Descendants
#' @importFrom phangorn Children
#' @importFrom combinat combn
#' @importFrom compiler cmpfun
#'
#' @examples
#'
#' ## generate a random tree
#' x <- rtree(6)
#'
#' ## create matrix of MRCAs: entry (i,j) is the node number of the MRCA of tips i and j
#' linear.mrca(x,6)
#'
linear.mrca <- function(tree,k=0) { # k is number of tips, which can be passed to the function to save on computation
  if (k==0) {k <- length(tree$tip.label)}
  M <- matrix(0, nrow=k, ncol=k); # initialise matrix
  T <- tree$Nnode # total number of internal nodes
  # traverse internal nodes from root down 
  for (tmp in (k+1):(k+T)){
    # find the children of tmp. Then tmp is the MRCA of all pairs of tips descending from different children
    tmp.desc <- Children(tree,tmp)
    Desc <- sapply(1:length(tmp.desc), function(x) Descendants(tree,tmp.desc[[x]],type="tips"))
      if (length(tmp.desc)==2) {  # tmp is the MRCA of tips descending from child one and tips from child two
        I <- Desc[[1]]; J <- Desc[[2]]
          for (i in I)  {
           for (j in J)  {
            M[i,j] <- M[j,i] <- tmp  
           }
          }
      }
      else { # for each pair of children of tmp, tmp is the MRCA of their descendant tips
        pairs <- combn(length(Desc),2)
        for (p in 1:length(pairs[1,])) {
          for (i in Desc[[pairs[1,p]]])  {
           for (j in Desc[[pairs[2,p]]])  {
            M[i,j] <- M[j,i] <- tmp  
           }
          }
        } 
      }
  }
  diag(M) <- 1:k # we define the diagonal elements of M to be the tips themselves
  return(M)
}
linear.mrca <- cmpfun(linear.mrca) # compile




#' Pendant edges
#'
#' Extract just the pendant edges from the vector \code{tree$edge}.
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree an object of the class \code{phylo} 
#' @param k number of tips in tree
#'
#' @importFrom compiler cmpfun
#'
pen.edge.tree <- function(tree,k) {tree$edge[match(1:k, tree$edge[,2]),] }
pen.edge.tree <- cmpfun(pen.edge.tree)



#' Pendant edges, matched
#'
#' Extract the pendant edges from the vector \code{tree$edge}, in the order given by \code{labelmatch}.
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree an object of the class \code{phylo} 
#' @param labelmatch a vector specifying the order of the tips in the output. This is used by other functions in the package to match the tip labels in the trees
#'
#' @importFrom compiler cmpfun
#'
pen.edge.treematch  <- function(tree,labelmatch) {tree$edge[match(labelmatch, tree$edge[,2]),] }
pen.edge.treematch <- cmpfun(pen.edge.treematch)





#' Tree vector function
#'
#' Function which takes a phylo as input and outputs the vector for the metric. 
#' The elements of the vector are numeric if \code{type="number"}, 
#' and otherwise they are functions of lambda.
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tr1 an object of the class \code{phylo} 
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{type="function"}. 
#' @param type logical which takes the inputs "\code{number}" (default) or "\code{function}". When \code{type="number"}, the entries of the output vector are numeric; when \code{type="function"} each entry is a function of lambda.
#'
#' @importFrom combinat combn2
#' @importFrom compiler cmpfun
#' @import ape 
#'
#' @examples
#'
#' ## generate a random tree
#' tree <- rtree(6)
#' ## topological vector of mrca distances from root:
#' tree.vec(tree) 
#' ## vector of mrca distances from root when lambda=0.5:
#' tree.vec(tree,0.5) 
#' ## vector of mrca distances as a function of lambda:
#' vec.func <- tree.vec(tree,type="function") 
#' ## evaluate the vector at lambda=0.5:
#' vec.func(0.5)  
#'
tree.vec <- function(tr1,lambda=0,type="number") { # allow output type to be number or function
  if (type=="number"){
    if (lambda<0) {stop("Pick lambda in [0,1]")}
    if (lambda>1) {stop("Pick lambda in [0,1]")}
    k <- length(tr1$tip.label)
    # checks and warnings
    
    if (lambda!=0) { # if lambda=0 then we don't need edge lengths to be defined, but if lambda!=0 then we do
      if (is.null(tr1$edge.length)) {
        stop("edge lengths not defined")
      }
    }
    
    M1 <- linear.mrca(tr1,k); # kxk MRCA matrix for tr1
    pairs <- combn2(1:k)
    tiporder <- order(tr1$tip.label)
        
    if (lambda!=1){ # make a copy with edge lengths = 1 because we need to know topological distances
      TR1 <- tr1;       TR1$edge.length <- rep(1,length(tr1$edge.length))
      D1 <- dist.nodes(TR1); 
    }
    if (lambda!=0) { # if lambda!=0 we need to know branch length distances
      d1 <- dist.nodes(tr1);
    }
    
    # vt is the purely topological vector (don't waste time computing if lambda=1)
    if (lambda==1) { vt <- rep(0,k*(k-1)/2)}
    else {
      vt <- apply(pairs, 1, function(x) D1[k+1,M1[[tiporder[[x[1]]],tiporder[[x[2]]]]]]) 
    }
    # append k entries of "1" for pendant edges
    vt <- as.numeric(c(vt,rep(1,k)))
    
    # vl is the purely length-based vector (don't waste time computing if lambda=0)
    if (lambda==0) { vl <- rep(0,k*(k+1)/2) }
    else {
      vl <- apply(pairs, 1, function(x) d1[k+1,M1[[tiporder[[x[1]]],tiporder[[x[2]]]]]]) 
      ep1 <- pen.edge.treematch(tr1,tiporder);
      pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])
      vl <- as.numeric(c(vl,pen.length1)) 
      }
    
    v <- (1-lambda)*vt + lambda*vl
  
  return(v)
  }
  if (type=="function") {  
    lambda <- integer()
    k <- length(tr1$tip.label)
    # checks and warnings
    if (is.null(tr1$edge.length)) {
      stop("edge lengths not defined")
    }
    
    M1 <- linear.mrca(tr1,k); # kxk MRCA matrix for tree 1
    pairs <- combn2(1:k)
    tiporder <- order(tr1$tip.label)
    
    # make a copy of the tree called TR1 with edge lengths = 1
    TR1 <- tr1
    TR1$edge.length <- rep(1,length(tr1$edge.length));
    D1 <- dist.nodes(TR1); 
    # find distances based on branch lengths:
    d1 <- dist.nodes(tr1);
    
    # vt is the purely topological vector, vl is the purely length-based vector 
    vt <- apply(pairs, 1, function(x) D1[k+1,M1[[tiporder[[x[1]]],tiporder[[x[2]]]]]]) 
    vl <- apply(pairs, 1, function(x) d1[k+1,M1[[tiporder[[x[1]]],tiporder[[x[2]]]]]]) 
    
    # append vector of pendant branch lengths
    ep1 <- pen.edge.treematch(tr1,tiporder);
    pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])
    
    vlambda <- function(lambda) {
      if (lambda<0) {stop("Pick lambda in [0,1]")}
      if (lambda>1) {stop("Pick lambda in [0,1]")}
      (c(((1-lambda)*vt + lambda*vl),(lambda*pen.length1))) }
    
    return(vlambda)
  }
}
tree.vec <- cmpfun(tree.vec)




#' Metric function
#'
#' Comparison of two trees using the Kendall Colijn metric
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tr1 an object of the class \code{phylo} 
#' @param tr2 an object of the class \code{phylo} 
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{type="function"}.
#' @param type logical which takes the inputs "\code{number}" (default) or "\code{function}". When \code{type="number"}, the output is numeric; when \code{type="function"} the output is a function of lambda.
#'
#'
#' @importFrom compiler cmpfun
#' @importFrom combinat combn2
#' @import ape 
#'
#'
#' @examples
#'
#' ## generate random trees
#' tree1 <- rtree(6); tree2 <- rtree(6)
#' tree.dist(tree1,tree2) # lambda=0
#' tree.dist(tree1,tree2,1)  # lambda=1
#' dist.func <- tree.dist(tree1,tree2,type="function") # distance as a function of lambda
#' dist.func(0) # evaluate at lambda=0. Equivalent to tree.dist(tree1,tree2).
#' ## We can see how the distance changes when moving from focusing on topology to length:
#' plot(sapply(seq(0,1,length.out=100), function(x) dist.func(x)), type="l",ylab="",xlab="") 
#'
#'
tree.dist <- function(tr1,tr2,lambda=0,type="number") { # allow output type to be number or function of lambda
  if (type=="number"){
    if (lambda<0) {stop("Pick lambda in [0,1]")}
    if (lambda>1) {stop("Pick lambda in [0,1]")}
    k <- length(tr1$tip.label)
    # checks and warnings
    if (k != length(tr2$tip.label)) {
      stop("trees have different numbers of tips")
    }
    if (setequal(tr1$tip.label,tr2$tip.label) == FALSE) {
      stop("trees have different tip label sets")
    } 
    if (lambda!=0) { # if lambda=0 then we don't need edge lengths to be defined, but if lambda!=0 then we do
      if (is.null(tr1$edge.length)) {
        stop("edge lengths not defined in first tree")
      }
      if (is.null(tr2$edge.length)) {
        stop("edge lengths not defined in second tree")
      }
    }
    
    M1 <- linear.mrca(tr1,k); # kxk MRCA matrix for tree 1
    M2 <- linear.mrca(tr2,k);
    labelmatch <- match(tr1$tip.label, tr2$tip.label);
    if (lambda!=1){ # make a copy of the trees called TR1 and TR2, with edge lengths = 1
      TR1 <- tr1; TR2 <- tr2
      TR1$edge.length <- rep(1,length(tr1$edge.length));
      TR2$edge.length <- rep(1,length(tr2$edge.length));
      D1 <- dist.nodes(TR1); # if lambda!=1 we need to know edge count distances
      D2 <- dist.nodes(TR2);
    }
    if (lambda!=0) { # if lambda!=0 we need to know branch length distances.
      d1 <- dist.nodes(tr1);
      d2 <- dist.nodes(tr2);
    }
    
    pairs <- combn2(1:k)
    # vt is the purely topological vector (don't waste time computing if lambda=1)
    # vl is the purely length-based vector (don't waste time computing if lambda=0)
    if (lambda==1) { vt <- rep(0,k*(k-1)/2)}
    else {
      vt <- apply(pairs, 1, function(x) D1[k+1,M1[[x[1],x[2]]]] - D2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]]) 
    }
    if (lambda==0) { vl <- rep(0,k*(k-1)/2)}
    else {
      vl <- apply(pairs, 1, function(x) d1[k+1,M1[[x[1],x[2]]]] - d2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]]) 
    }
    
    v <- (1-lambda)*vt + lambda*vl
    
    if (lambda!=0) {
      # append vector of difference in pendant branch lengths
      ep1 <- pen.edge.tree(tr1,k);
      ep2 <- pen.edge.treematch(tr2,labelmatch);
      pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])
      pen.length2 <- apply(ep2, 1, function(x) d2[x[1],x[2]])
      pen.length.diff <- sapply(1:k, function(x) pen.length1[[x]] - pen.length2[[x]])
      v <- as.numeric(c(v,lambda*pen.length.diff)) 
    }
    
    return(sqrt(sum(v^2)))
  }
  if (type=="function") {  
    lambda <- integer()
    k <- length(tr1$tip.label)
    # checks and warnings
    if (k != length(tr2$tip.label)) {
      stop("trees have different numbers of tips")
    }
    if (setequal(tr1$tip.label,tr2$tip.label) == FALSE) {
      stop("trees have different tip label sets")
    } 
    if (is.null(tr1$edge.length)) {
      stop("edge lengths not defined in first tree")
    }
    if (is.null(tr2$edge.length)) {
      stop("edge lengths not defined in second tree")
    }
    
    M1 <- linear.mrca(tr1,k); # kxk MRCA matrix for tree 1
    M2 <- linear.mrca(tr2,k);
    labelmatch <- match(tr1$tip.label, tr2$tip.label);
    # make a copy of the trees called TR1 and TR2, with edge lengths = 1
    TR1 <- tr1; TR2 <- tr2
    TR1$edge.length <- rep(1,length(tr1$edge.length));
    TR2$edge.length <- rep(1,length(tr2$edge.length));
    D1 <- dist.nodes(TR1); 
    D2 <- dist.nodes(TR2);
    
    # get full distance matrices with lengths
    d1 <- dist.nodes(tr1);
    d2 <- dist.nodes(tr2);
    
    pairs <- combn2(1:k)
    # vt is the purely topological vector, vl is the purely length-based vector 
    vt <- apply(pairs, 1, function(x) D1[k+1,M1[[x[1],x[2]]]] - D2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]]) 
    vl <- apply(pairs, 1, function(x) d1[k+1,M1[[x[1],x[2]]]] - d2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]]) 
    
    # append vector of difference in pendant branch lengths
    ep1 <- pen.edge.tree(tr1,k);
    ep2 <- pen.edge.treematch(tr2,labelmatch);
    pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])
    pen.length2 <- apply(ep2, 1, function(x) d2[x[1],x[2]])
    pen.length.diff <- sapply(1:k, function(x) pen.length1[[x]] - pen.length2[[x]])
    
    vlambda <- function(lambda) {
      if (lambda<0) {stop("Pick lambda in [0,1]")}
      if (lambda>1) {stop("Pick lambda in [0,1]")}
      sqrt(sum((c(((1-lambda)*vt + lambda*vl),(lambda*pen.length.diff)))^2)) }
    
    return(vlambda)
  }
}
tree.dist <- cmpfun(tree.dist)



#' Metric function for \code{multiPhylo} input
#'
#' Comparison of a list of trees using the Kendall Colijn metric. Output is given as a pairwise distance matrix.
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param trees an object of the class \code{multiPhylo} 
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{type="function"}. 
#' @param type logical which takes the inputs "\code{number}" (default) or "\code{function}". When \code{type="number"}, the output is a matrix where each entry is numeric; when \code{type="function"} the output is a matrix where each entry is a function of lambda.
#'
#'
#' @importFrom compiler cmpfun
#' @importFrom combinat combn2
#' @import ape 
#'
#'
#' @examples
#'
#' ## generate 10 random trees, each with 6 tips
#' trees <- rmtree(10,6)
#' 
#' ## pairwise distance matrix when lambda=0
#' multi.dist(trees)
#' 
#' ## pairwise distance matrix as a function of lambda:
#' m <- multi.dist(trees, type="function")
#' 
#' ## evaluate at lambda=0. Equivalent to multi.dist(trees).
#' m(0)
#'
#' ## A method to visualise these distances with MDS:
#' require(ade4)
#' 
#' ## find an optimum projection of the points in 2 dimensions:
#' mMDS <- dudi.pco(as.dist(m(0)), scannf=FALSE,nf=2)
#' 
#' ## put the coordinates of these points into a data frame
#' mdf <-mMDS$li
#'
#' ## basic ade4 plot
#' s.label(mdf)
#' 
#' ## ggplot2 version
#' if(require(ggplot2) && require(RColorBrewer)){
#' mplot <- ggplot(mdf, aes(x=A1, y=A2)) 
#' 
#' mpalette <- brewer.pal(10,"Paired") # create colour palette
#' 
#' mplot + geom_point(colour=mpalette,size=5) +  
#'   xlab("") + ylab("") + theme_bw(base_family = "")
#' }
#' 
#' 
#' ## An example using data:
#' ## These woodmice phylogenies were created using the bootstrapping example in package \code{ape}
#' data(woodmiceTrees)
#' woodmiceDists <- multi.dist(woodmiceTrees) # find topological distances
#' woodmiceMDS <- dudi.pco(as.dist(woodmiceDists), scannf=FALSE, nf=2)
#' woodmicedf <- woodmiceMDS$li
#'
#' if(require(ggplot2)){
#' woodmiceplot <- ggplot(woodmicedf, aes(x=A1, y=A2)) # create plot
#' woodmiceplot + geom_density2d(colour="gray80") + # contour lines
#'  geom_point(size=6, shape=1, colour="gray50") + # grey edges
#'  geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
#'  xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
#' }
#'
#' \dontrun{
#' if(require(rgl)){
#' woodmiceMDS3D <- dudi.pco(as.dist(woodmiceDists), scannf=FALSE, nf=3)
#' plot3d(woodmiceMDS3D$li[,1], woodmiceMDS3D$li[,2], woodmiceMDS3D$li[,3], type="s", size=1.5, 
#'    col="navy", alpha=0.5, xlab="", ylab="", zlab="")
#' }
#' }
multi.dist <- function(trees,lambda=0,type="number") { # allow output type to be number or function
  #checks and warnings
  if (class(trees) != "multiPhylo"){
    stop("input must be of class multiPhylo")
  }
  l <- length(trees)
  k <- length(trees[[1]]$tip.label)
  
  for (i in 1:l) {
    if (k != length(trees[[i]]$tip.label)) {
      stop("trees must all have the same number of tips")
    }
    if (setequal(trees[[i]]$tip.label,trees[[1]]$tip.label) == FALSE) {
      stop("trees have different tip label sets")
    } 
  }
  if (type=="number"){ 
    # checks and warnings  
    if (lambda<0) {stop("Pick lambda in [0,1]")}
    if (lambda>1) {stop("Pick lambda in [0,1]")}
    if (lambda!=0) { # if lambda=0 then we don't need edge lengths to be defined, but if lambda!=0 then we do
      if (is.null(trees[[i]]$edge.length)) {
        stop("edge lengths not defined")
      }}
    
    M <- lapply(1:l, function(x) linear.mrca(trees[[x]],k));
    labelmatch <- lapply(1:l, function (y) 
      match(trees[[1]]$tip.label,trees[[y]]$tip.label));
    pairs <- combn2(1:k)
    x <- k*(k-1)/2
    
    # topvecs is the purely topological matrix of vectors (don't waste time computing if lambda=1)
    # lvecs is the purely length-based matrix of vectors (don't waste time computing if lambda=0)
    if (lambda!=1) {
      # make a copy of the trees with edge lengths = 1
      TREES <- trees
      for (i in 1:l) {
        TREES[[i]]$edge.length <- rep(1,length(trees[[i]]$edge.length));
      }
      D <- sapply(1:l, function(x) dist.nodes(TREES[[x]])[k+1,]); # vector of vectors
      
      topvecs <- sapply(1:l, function(y) apply(pairs, 1, function(x) D[M[[y]][[labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]],y])); 
      tv <- (1-lambda)*topvecs
    }
    if (lambda!=0) {
      # we also need to know branch length distance matrix
      d <- lapply(1:l, function(x) dist.nodes(trees[[x]]));
      
      lvecs <- sapply(1:l, function(y) apply(pairs, 1, function(x) d[[y]][k+1,M[[y]][labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]])); 
      tl <- lambda*lvecs  
      
      # append vector of difference in pendant branch lengths
      E <- lapply(1:l, function(x) pen.edge.treematch(trees[[x]],labelmatch[[x]]))
      Pen <- sapply(1:l, function(x) apply(E[[x]], 1, function(y) d[[x]][y[1],y[2]]) )
      P <- lambda*Pen
    }
    
    if (lambda==0) {
      # matrix where each entry is a vector of squared differences
      sqdistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
        if (a<b) sapply(1:x, function(n) (tv[[n,a]]-tv[[n,b]])^2)))
    }
    else if (lambda==1) {
      # matrix where each entry is a vector of squared differences plus pendant lengths
      sqdistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
        if (a<b) c(sapply(1:x, function(n) (tl[[n,a]]-tl[[n,b]])^2),
                   sapply(1:k, function(d) (P[,a][[d]]-P[,b][[d]])^2))))
    }
    else {   
      # matrix where each entry is a vector of squared differences plus pendant lengths
      sqdistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
        if (a<b) c(sapply(1:x, function(n) (tv[[n,a]]+tl[[n,a]]-tv[[n,b]]-tl[[n,b]])^2),
                   sapply(1:k, function(d) (P[,a][[d]]-P[,b][[d]])^2))))
    }
    # final matrix: each entry is square root of sum of sqdistmat entry
    distmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
      sqrt(sum(sqdistmat[[b,a]]))))
    return(as.dist(distmat))
  } 
  else if (type=="function"){ 
    lambda <- integer()
    # check: we need edge lengths defined
    if (is.null(trees[[i]]$edge.length)) {
      stop("edge lengths not defined")
    }
    
    
    M <- lapply(1:l, function(x) linear.mrca(trees[[x]],k));
    labelmatch <- lapply(1:l, function (y) 
      match(trees[[1]]$tip.label,trees[[y]]$tip.label));
    pairs <- combn2(1:k)
    x <- k*(k-1)/2
    
    # make a copy of the trees with edge lengths = 1
    TREES <- trees
    for (i in 1:l) {
      TREES[[i]]$edge.length <- rep(1,length(trees[[i]]$edge.length));
    }
    D <- sapply(1:l, function(x) dist.nodes(TREES[[x]])[k+1,]); # vector of vectors
    
    tv <- sapply(1:l, function(y) apply(pairs, 1, function(x) D[M[[y]][[labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]],y])); 
    
    # we also need to know branch length distance matrix
    d <- lapply(1:l, function(x) dist.nodes(trees[[x]]));
    
    tl <- sapply(1:l, function(y) apply(pairs, 1, function(x) d[[y]][k+1,M[[y]][labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]])); 
    
    # append vector of difference in pendant branch lengths
    E <- lapply(1:l, function(x) pen.edge.treematch(trees[[x]],labelmatch[[x]]))
    P <- sapply(1:l, function(x) apply(E[[x]], 1, function(y) d[[x]][y[1],y[2]]) )
    
    # matrix where each entry is a vector of squared differences plus pendant lengths
    sqdistmat <- function(lambda) {as.dist(sapply(1:l, function(a) sapply(1:l, function(b) 
      if (a<b) {sqrt(sum(c(sapply(1:x, function(n) ((1-lambda)*(tv[[n,a]]-tv[[n,b]])+lambda*(tl[[n,a]]-tl[[n,b]]))^2),
                           sapply(1:k, function(d) (lambda*(P[,a][[d]]-P[,b][[d]]))^2))))}
      else {0})))  }
    
  } 
  return(sqdistmat)
}
multi.dist <- cmpfun(multi.dist)


#' Geometric median tree function
#'
#' Finds the geometric median of a set of trees according to the Kendall Colijn metric.
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param trees an object of the class \code{multiPhylo} 
#' @param likes a vector of weightings for the trees. Defaults to a vector of 1's so that all trees are equally weighted, but can be used to weight trees according to likelihood or other characteristics.
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{type="function"}.
#' 
#'
#' @importFrom compiler cmpfun
#' @importFrom combinat combn2
#' @import ape 
#'
#'
#' @examples
#'
#' ## generate 10 random trees, each with 6 tips
#' trees <- rmtree(10,6)
#' ## Geometric median tree:
#' mymedian <- med.tree(trees)
#' mymedian$centre # the vector at the 'centre' of the trees; may not correspond to an actual tree
#' mymedian$median # the identifier(s) of the tree(s) closest to the central vector
#' mymedian$mindist # the distance of the median tree(s) from the central vector
#'
#' ## Example with woodmice data:
#' data(woodmiceTrees)
#' woodmiceMed <- med.tree(woodmiceTrees)
#' ## plot the (first) geometric median tree (there are seven topologically identical median trees):
#' plot(woodmiceTrees[[woodmiceMed$median[[1]]]])
#' 
#' ## finding the geometric median tree from a single cluster:
#' woodmiceDists <- multi.dist(woodmiceTrees)
#' woodmiceMDS <- dudi.pco(as.dist(woodmiceDists), scannf=FALSE, nf=2)
#' ## isolate the trees from the largest cluster
#' woodmiceCluster1 <- woodmiceTrees[intersect(
#'    intersect(which(woodmiceMDS$li[,1]>(-2)),which(woodmiceMDS$li[,1]<2)),
#'    intersect(which(woodmiceMDS$li[,2]>(-2.5)),which(woodmiceMDS$li[,1]<2.5))
#'    )]
#' ## find the geometric median
#' geomMedWoodmice1 <- med.tree(woodmiceCluster1)
#' plot(woodmiceCluster1[[geomMedWoodmice1$median[[1]]]])
#' # this has the same topology as the overall median tree:
#' tree.dist(woodmiceTrees[[woodmiceMed$median[[1]]]],
#'   woodmiceCluster1[[geomMedWoodmice1$median[[1]]]])
#'
#' ## However, median trees from other clusters have different topologies, for example:
#' ## isolate the trees from the second largest cluster:
#' woodmiceCluster2 <- woodmiceTrees[intersect(
#'  intersect(which(woodmiceMDS$li[,1]>(-1)),which(woodmiceMDS$li[,1]<8)),
#'  intersect(which(woodmiceMDS$li[,2]>1),which(woodmiceMDS$li[,1]<6))
#' )]
#' ## find the geometric median
#' geomMedWoodmice2 <- med.tree(woodmiceCluster2)
#' plot(woodmiceCluster2[[geomMedWoodmice2$median[[1]]]])
#' ## This is another representative topology, which is different from those we found above:
#' tree.dist(woodmiceCluster2[[geomMedWoodmice2$median[[1]]]],
#'   woodmiceCluster2[[geomMedWoodmice2$median[[1]]]])
#'
med.tree <- function(trees,likes=rep(1,length(trees)),lambda=0) {
  n <- length(trees)
  if (length(likes)!=n) {stop("Number of likelihoods is not equal to number of trees.")}
  if (lambda<0) {stop("Pick lambda in [0,1]")}
  if (lambda>1) {stop("Pick lambda in [0,1]")}
  
  k <- length(trees[[1]]$tip.label)
  for (i in 1:n) {
    if (k != length(trees[[i]]$tip.label)) {
      stop("trees must all have the same number of tips")
    }
    if (setequal(trees[[i]]$tip.label,trees[[1]]$tip.label) == FALSE) {
      stop("trees have different tip label sets")
    } 
  }
  if (lambda!=0) { # if lambda=0 then we don't need edge lengths to be defined, but if lambda!=0 then we do
    if (is.null(trees[[i]]$edge.length)) {
      stop("edge lengths not defined")
    }}
  
  labelmatch <- lapply(1:n, function (y) 
    match(trees[[1]]$tip.label,trees[[y]]$tip.label))
  
  # version of tree.vec which applies labelmatch first
  tree.vec.match <- function(tr1,lambda,labelmatchi,k,n) { 
    M1 <- linear.mrca(tr1,k); # kxk MRCA matrix for tree 
    
    if (lambda!=1){ # make a copy with edge lengths = 1
      TR1 <- tr1
      TR1$edge.length <- rep(1,length(tr1$edge.length))
      D1 <- dist.nodes(TR1)
    }
    if (lambda!=0) { # if lambda!=0 we need to know branch length distances
      # first, we need to rescale branch lengths so median is 1
      tr1$edge.length <- tr1$edge.length/median(tr1$edge.length)
      d1 <- dist.nodes(tr1);
    }
    
    pairs <- combn2(1:k)
    # vt is the purely topological vector (don't waste time computing if lambda=1)
    # vl is the purely length-based vector (don't waste time computing if lambda=0)
    if (lambda==1) { vt <- rep(0,k*(k-1)/2)}
    else {
      vt <- apply(pairs, 1, function(x) D1[k+1,M1[[labelmatchi[x[1]],labelmatchi[x[2]]]]]) 
    }
    if (lambda==0) { vl <- rep(0,k*(k-1)/2)}
    else {
      vl <- apply(pairs, 1, function(x) d1[k+1,M1[[labelmatchi[x[1]],labelmatchi[x[2]]]]]) 
    }
    
    v <- (1-lambda)*vt + lambda*vl
    
    if (lambda!=0) {
      # append vector of pendant branch lengths
      ep1 <- pen.edge.treematch(tr1,labelmatchi)
      pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])
      v <- as.numeric(c(v,lambda*pen.length1)) 
    }
    
    return(v)
  }
  # initialise vector, length n choose 2 if lambda=0, otherwise n choose 2 + n
  if (lambda==0) {  centre <- rep(0,(k*(k-1)/2)) }
  else {  centre <- rep(0,(k*(k+1)/2)) }
  vecs <- list()
  for (i in 1:n) {
    vecs[[i]] <- tree.vec.match(trees[[i]],lambda,labelmatch[[i]],k,n)
    centre <- centre + vecs[[i]]*likes[[i]]
  }
  centre <- centre/n
  # also want to know which vecs[[i]] is closest to median
  d <- list()
  for (i in 1:n){
    v <- vecs[[i]]-centre
    d[[i]] <- sqrt(sum(v^2))
  }
  class(d) <- "numeric"
  md <- min(d)
  median <- which(d==md)
  
  result <- list()
  result$centre <- centre
  result$median <- median
  result$mindist <- md
  
  return(result)
}
med.tree <- cmpfun(med.tree)
