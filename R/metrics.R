
#' Function to make mrca matrix of a tree, where entry (i,j) gives the mrca of tips i and j
#'
#' Description of this function..
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @export
#'
#' @param tree ...
#' @param k ...

#' @import ape
#' @import phangorn
#' @import compiler
#' @import fastmatch
#' @import combinat
#'
linear.mrca <- function(tree,k)
{
    M <- matrix(0, nrow=k, ncol=k); # initialise matrix
                                        # traverse internal nodes from root down
    for (tmp in (k+1):(2*k-1)){
                                        # find the two children of tmp
        tmp.desc <- Children(tree,tmp)
                                        # tmp is the MRCA of all pairs of tips descending from child one and child two
        I <- Descendants(tree,tmp.desc[[1]], type="tips")
        J <- Descendants(tree,tmp.desc[[2]], type="tips")
        for (i in I)  {
            for (j in J)  {
                M[i,j] <- M[j,i] <- tmp
            }  }
    }
    return(M)
}
linear.mrca <- cmpfun(linear.mrca) # compile




#' Title of the function
#'
#' Description of this function..
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree ...
#' @param k ...
#'
#'
pen.edge.tree <- function(tree,k) {tree$edge[fmatch(1:k, tree$edge[,2]),] }
pen.edge.tree <- cmpfun(pen.edge.tree)



#' Title of the function
#'
#' Description of this function..
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree ...
#' @param labelmatch ...
#'
pen.edge.treematch  <- function(tree,labelmatch) {tree$edge[fmatch(labelmatch, tree$edge[,2]),] }
pen.edge.treematch <- cmpfun(pen.edge.treematch)





#' Title of the function
#'
#' Description of this function..
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param x an object of the class \code{phylo} ...
#' @param lambda ...
#' @param type ...
#'
CK.metric <- function(x,lambda=0,type="number") { # allow output type to be number or function
    if (type=="number"){
        ## checks and warnings
        if (lambda<0) {stop("Pick lambda in [0,1]")}
        if (lambda>1) {stop("Pick lambda in [0,1]")}
        k <- length(x$tip.label)

        if (lambda!=0) { # if lambda=0 then we don't need edge lengths to be defined, but if lambda!=0 then we do
            if (is.null(x$edge.length)) {
                stop("edge lengths not defined")
            }
        }

        M1 <- linear.mrca(x,k); # kxk MRCA matrix for tree 1

        if (lambda!=1){ # make a copy with edge lengths = 1
            X <- x
            X$edge.length <- rep(1,2*k-2);
            D1 <- dist.nodes(X); # if lambda!=1 we need to know edge count distances
        }
        if (lambda!=0) { # if lambda!=0 we need to know branch length distances
            d1 <- dist.nodes(x);
        }

        pairs <- combn2(1:k)
                                        # vt is the purely topological vector (don't waste time computing if lambda=1)
                                        # vl is the purely length-based vector (don't waste time computing if lambda=0)
        if (lambda==1) { vt <- rep(0,k*(k-1)/2)}
        else {
            vt <- apply(pairs, 1, function(x) D1[k+1,M1[[x[1],x[2]]]])
        }
        if (lambda==0) { vl <- rep(0,k*(k-1)/2)}
        else {
            vl <- apply(pairs, 1, function(x) d1[k+1,M1[[x[1],x[2]]]])
        }

        v <- (1-lambda)*vt + lambda*vl

        if (lambda!=0) {
                                        # append vector of pendant branch lengths
            ep1 <- pen.edge.tree(x,k);
            pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])
            v <- as.numeric(c(v,lambda*pen.length1))
        }

        return(v)
    }
    if (type=="function") {
        lambda <- integer()
        k <- length(x$tip.label)
                                        # checks and warnings
        if (is.null(x$edge.length)) {
            stop("edge lengths not defined")
        }

        M1 <- linear.mrca(x,k); # kxk MRCA matrix for tree 1

        ## make a copy of the tree called X with edge lengths = 1
        X <- x
        X$edge.length <- rep(1,2*k-2);
        D1 <- dist.nodes(X);
        ## find distances based on branch lengths:
        d1 <- dist.nodes(x);

        pairs <- combn2(1:k)
        ## vt is the purely topological vector, vl is the purely length-based vector
        vt <- apply(pairs, 1, function(x) D1[k+1,M1[[x[1],x[2]]]])
        vl <- apply(pairs, 1, function(x) d1[k+1,M1[[x[1],x[2]]]])

        ## append vector of pendant branch lengths
        ep1 <- pen.edge.tree(x,k);
        pen.length1 <- apply(ep1, 1, function(x) d1[x[1],x[2]])

        vlambda <- function(lambda) {
            if (lambda<0) {stop("Pick lambda in [0,1]")}
            if (lambda>1) {stop("Pick lambda in [0,1]")}
            (c(((1-lambda)*vt + lambda*vl),(lambda*pen.length1))) }

        return(vlambda)
    }
}
CK.metric <- cmpfun(CK.metric)




#' Title of the function
#'
#' Description of this function..
#'
#' @export
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tr1 ...
#' @param tr2 ...
#' @param lambda ...
#' @param type ...
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
        labelmatch <- fmatch(tr1$tip.label, tr2$tip.label);
        if (lambda!=1){ # make a copy of the trees called TR1 and TR2, with edge lengths = 1
            TR1 <- tr1; TR2 <- tr2
            TR1$edge.length <- rep(1,2*k-2);
            TR2$edge.length <- rep(1,2*k-2);
            D1 <- dist.nodes(TR1); # if lambda!=1 we need to know edge count distances
            D2 <- dist.nodes(TR2);
        }
        if (lambda!=0) { # if lambda!=0 we need to know branch length distances.
            d1 <- dist.nodes(tr1);
            d2 <- dist.nodes(tr2);
        }

        pairs <- combn2(1:k)
        ## vt is the purely topological vector (don't waste time computing if lambda=1)
        ## vl is the purely length-based vector (don't waste time computing if lambda=0)
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
            ## append vector of difference in pendant branch lengths
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
        labelmatch <- fmatch(tr1$tip.label, tr2$tip.label);
        ## make a copy of the trees called TR1 and TR2, with edge lengths = 1
        TR1 <- tr1; TR2 <- tr2
        TR1$edge.length <- rep(1,2*k-2);
        TR2$edge.length <- rep(1,2*k-2);
        D1 <- dist.nodes(TR1);
        D2 <- dist.nodes(TR2);

        ## get full distance matrices with lengths
        d1 <- dist.nodes(tr1);
        d2 <- dist.nodes(tr2);

        pairs <- combn2(1:k)
        ## vt is the purely topological vector, vl is the purely length-based vector
        vt <- apply(pairs, 1, function(x) D1[k+1,M1[[x[1],x[2]]]] - D2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]])
        vl <- apply(pairs, 1, function(x) d1[k+1,M1[[x[1],x[2]]]] - d2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]])

        ## append vector of difference in pendant branch lengths
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






## multi.dist <- function(trees,lambda=0,type="number") { # allow output type to be number or function
##     ##checks and warnings
##     if (class(trees) != "multiPhylo"){
##         stop("input must be of class multiPhylo")
##     }
##     l <- length(trees)
##     k <- length(trees[[1]]$tip.label)

##     for (i in 1:l) {
##         if (k != length(trees[[i]]$tip.label)) {
##             stop("trees must all have the same number of tips")
##         }
##         if (setequal(trees[[i]]$tip.label,trees[[1]]$tip.label) == FALSE) {
##             stop("trees have different tip label sets")
##         }
##     }
##     if (type=="number"){
##         ## checks and warnings
##         if (lambda<0) {stop("Pick lambda in [0,1]")}
##         if (lambda>1) {stop("Pick lambda in [0,1]")}
##         if (lambda!=0) { # if lambda=0 then we don't need edge lengths to be defined, but if lambda!=0 then we do
##             if (is.null(trees[[i]]$edge.length)) {
##                 stop("edge lengths not defined")
##             }}

##         M <- lapply(1:l, function(x) linear.mrca(trees[[x]],k));
##         labelmatch <- lapply(1:l, function (y)
##                              fmatch(trees[[1]]$tip.label,trees[[y]]$tip.label));
##         pairs <- combn2(1:k)
##         x <- k*(k-1)/2

##         ## topvecs is the purely topological matrix of vectors (don't waste time computing if lambda=1)
##         ## lvecs is the purely length-based matrix of vectors (don't waste time computing if lambda=0)
##         if (lambda!=1) {
##             ## make a copy of the trees with edge lengths = 1
##             TREES <- trees
##             for (i in 1:l) {
##                 TREES[[i]]$edge.length <- rep(1,(2*k-2));
##             }
##             D <- sapply(1:l, function(x) dist.nodes(TREES[[x]])[k+1,]); # vector of vectors

##             topvecs <- sapply(1:l, function(y) apply(pairs, 1, function(x) D[M[[y]][[labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]],y]));
##             tv <- (1-lambda)*topvecs
##         }
##         if (lambda!=0) {
##             ## we also need to know branch length distance matrix
##             d <- lapply(1:l, function(x) dist.nodes(trees[[x]]));

##             lvecs <- sapply(1:l, function(y) apply(pairs, 1, function(x) d[[y]][k+1,M[[y]][labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]]));
##             tl <- lambda*lvecs

##             ## append vector of difference in pendant branch lengths
##             E <- lapply(1:l, function(x) pen.edge.treematch(trees[[x]],labelmatch[[x]]))
##             Pen <- sapply(1:l, function(x) apply(E[[x]], 1, function(y) d[[x]][y[1],y[2]]) )
##             P <- lambda*Pen
##         }

##         if (lambda==0) {
##             ## matrix where each entry is a vector of squared differences
##             sqdistmat <- sapply(1:l, function(a) sapply(1:l, function(b)
##                                                         if (a<b) sapply(1:x, function(n) (tv[[n,a]]-tv[[n,b]])^2)))
##         }
##         else if (lambda==1) {
##             ## matrix where each entry is a vector of squared differences plus pendant lengths
##             sqdistmat <- sapply(1:l, function(a) sapply(1:l, function(b)
##                                                         if (a<b) c(sapply(1:x, function(n) (tl[[n,a]]-tl[[n,b]])^2),
##                                                                    sapply(1:k, function(d) (P[,a][[d]]-P[,b][[d]])^2))))
##         }
##         else {
##             ## matrix where each entry is a vector of squared differences plus pendant lengths
##             sqdistmat <- sapply(1:l, function(a) sapply(1:l, function(b)
##                                                         if (a<b) c(sapply(1:x, function(n) (tv[[n,a]]+tl[[n,a]]-tv[[n,b]]-tl[[n,b]])^2),
##                                                                    sapply(1:k, function(d) (P[,a][[d]]-P[,b][[d]])^2))))
##         }
##         ## final matrix: each entry is square root of sum of sqdistmat entry
##         distmat <- sapply(1:l, function(a) sapply(1:l, function(b)
##                                                   sqrt(sum(sqdistmat[[b,a]]))))
##         return(distmat)
##     }
##     else if (type=="function"){
##         lambda <- integer()
##         ## check: we need edge lengths defined
##         if (is.null(trees[[i]]$edge.length)) {
##             stop("edge lengths not defined")
##         }}


##     M <- lapply(1:l, function(x) linear.mrca(trees[[x]],k));
##     labelmatch <- lapply(1:l, function (y)
##                          fmatch(trees[[1]]$tip.label,trees[[y]]$tip.label));
##     pairs <- combn2(1:k)
##     x <- k*(k-1)/2

##     ## make a copy of the trees with edge lengths = 1
##     TREES <- trees
##     for (i in 1:l) {
##         TREES[[i]]$edge.length <- rep(1,(2*k-2));
##     }
##     D <- sapply(1:l, function(x) dist.nodes(TREES[[x]])[k+1,]); # vector of vectors

##     tv <- sapply(1:l, function(y) apply(pairs, 1, function(x) D[M[[y]][[labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]],y]));

##     ## we also need to know branch length distance matrix
##     d <- lapply(1:l, function(x) dist.nodes(trees[[x]]));

##     tl <- sapply(1:l, function(y) apply(pairs, 1, function(x) d[[y]][k+1,M[[y]][labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]]));

##     ## append vector of difference in pendant branch lengths
##     E <- lapply(1:l, function(x) pen.edge.treematch(trees[[x]],labelmatch[[x]]))
##     P <- sapply(1:l, function(x) apply(E[[x]], 1, function(y) d[[x]][y[1],y[2]]) )

##     ## matrix where each entry is a vector of squared differences plus pendant lengths
##     sqdistmat <- function(lambda) { sapply(1:l, function(a) sapply(1:l, function(b)
##                                                                    if (a<b) {sqrt(sum(c(sapply(1:x, function(n) ((1-lambda)*(tv[[n,a]]-tv[[n,b]])+lambda*(tl[[n,a]]-tl[[n,b]]))^2),
##                                                                                         sapply(1:k, function(d) (lambda*(P[,a][[d]]-P[,b][[d]]))^2))))}
##                                                                    else {0}))  }

## }
## return(sqdistmat)
## }
## multi.dist <- cmpfun(multi.dist)
