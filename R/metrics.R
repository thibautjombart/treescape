
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



#' Tree vector function
#'
#' Function which takes an object of class phylo as input and outputs the vector for the metric. 
#' The elements of the vector are numeric if \code{return_lambda_function=FALSE} (default), 
#' and otherwise they are functions of lambda.
#'
#' @export
#'
#' @author Jacob Almagro-Garcia \email{nativecoder@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree an object of the class \code{phylo}
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{type="function"}. 
#' @param return_lambda_function If true, a function that can be invoked with different lambda values is returned. This function returns the vector of metric values for the given lambda.
#'
#' @return The vector with the metric values or a function that produces the vector given a value of lambda.
#'
#' @import Rcpp
#' @import inline
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
#' vec.func <- tree.vec(tree,return_lambda_function=TRUE) 
#' ## evaluate the vector at lambda=0.5:
#' vec.func(0.5)  
#'
tree.vec <- function(tree, lambda=0, return_lambda_function=F) {

  if(lambda<0 || lambda>1) stop("Pick lambda in [0,1]")
  
  num_leaves <- length(tree$tip.label)
  num_edges <- nrow(tree$edge)
  
  # We work with ordered labels, using this vector to transform indices.
  tip_order <- match(1:num_leaves, order(tree$tip.label))
  
  # Ordering the edges by first column places the root at the bottom.
  # Descendants will be placed always before parents.
  edge_order <- order(tree$edge[,1], decreasing=T)
  edges <- tree$edge[edge_order,]
  edge_lengths <- tree$edge.length[edge_order]
  
  # We annotated the nodes of the tree in this list. In two passes we are going to
  # compute the partition each node induces in the tips (bottom-up pass) and the distance 
  # (in branch length and number of branches) from the root to each node (top-down pass).
  annotated_nodes <- list()
  
  # Bottom up (compute partitions, we store the branch lengths to compute distances 
  # to the root on the way down).
  for(i in 1:num_edges) {
      
    parent <- edges[i,1]
    child <- edges[i,2]
    
    # Initialization (leaves).
    if(child <= num_leaves) {
      # We translate the index for the sorted labels.
      child <- tip_order[child]
      # Leaves have as children themselves.
      annotated_nodes[[child]] <- list(root_distance=NULL, edges_to_root=1, partitions=list(child)) 
    }
    
    # Aggregate the children partitions (only if we are not visiting a leaf).
    aggregated_partitions <- annotated_nodes[[child]]$partitions[[1]]
    if((child > num_leaves)) {
      for(p in 2:length(annotated_nodes[[child]]$partitions))
        aggregated_partitions <- c(aggregated_partitions, annotated_nodes[[child]]$partitions[[p]])
    }
    
    # Update the branch length on the child.
    annotated_nodes[[child]]$root_distance <- edge_lengths[i]
    
    # We have not visited this internal node before.
    if(parent > length(annotated_nodes) || is.null(annotated_nodes[[parent]])) {
      # Assume the first time we get the left child partition.
      annotated_nodes[[parent]] <- list(root_distance=NULL, edges_to_root=1, partitions=list(aggregated_partitions))
    }
    # This is not the first time we have visited the node.
    else {
      # We store the next partition of leaves.
      annotated_nodes[[parent]]$partitions[[length(annotated_nodes[[parent]]$partitions)+1]] <- aggregated_partitions
    }
  }
  
  # Update the distance to the root at the root (i.e. 0)
  # And the number of edges to the root (i.e. 0).
  annotated_nodes[[num_leaves+1]]$root_distance <- 0
  annotated_nodes[[num_leaves+1]]$edges_to_root <- 0
  
  # Top down, compute distances to the root for each node.
  for(i in num_edges:1) {
    parent <- edges[i,1]
    child <- edges[i,2]
    
    # If the child is a leaf we translate the index for the sorted labels.
    if(child <= num_leaves)
      child <- tip_order[child]
    
    annotated_nodes[[child]]$root_distance <- annotated_nodes[[child]]$root_distance + annotated_nodes[[parent]]$root_distance
    annotated_nodes[[child]]$edges_to_root <- annotated_nodes[[child]]$edges_to_root + annotated_nodes[[parent]]$edges_to_root
  }
  
  # Distance vectors
  vector_length <- (num_leaves*(num_leaves-1)/2) + num_leaves
  length_root_distances <- double(vector_length)
  topological_root_distances <- integer(vector_length)
  
  # Fill-in the leaves (notice the involved index translation for leaves).
  topological_root_distances[(vector_length-num_leaves+1):vector_length] <- 1
  length_root_distances[(vector_length-num_leaves+1):vector_length] <- edge_lengths[match(1:num_leaves, edges[,2])][order(tree$tip.label)]
  
  # Instead of computing the lexicographic order for the combination pairs assume we
  # are filling in a symmetric distance matrix (using only the triangular upper part).
  # We just need to "roll" the matrix indices into the vector indices.
  # Examples for (k=5)
  # The combination c(1,4) would be located at position 3 on the vector.
  # The combination c(2,1) would be located at position 1 on the vector because d(2,1) = d(1,2).
  # The combination c(2,3) would be located at position 5 on the vector.
  
  index_offsets <- c(0, cumsum((num_leaves-1):1))
  
  # This is the slow part, we compute both vectors as gain would be marginal.
  sapply(annotated_nodes, function(node) {
    
    # We skip leaves and the root (if the MRCA for M groups of leaves is at the root
    # all combinations of leaves -among different groups- have 0 as distance to the root).
    # For large trees this can spare us of computing a lot of combinations. 
    # Example: In a perfectly balanced binary tree (N/2 leaves at each side of the root), 
    # at the root we'd save (N/2) * (N/2) combinations to update. Worst case scenario is completely 
    # unbalanced tree (N-1,1), we'd save in that case only N-1 combinations.
    
    # Make sure we are not visiting a leaf or the root.
    if(length(node$partitions) > 1 && node$root_distance > 0) {
      
      # Update all combinations for pairs of leaves from different groups.
      num_groups <- length(node$partitions)
      for(group_a in 1:(num_groups-1)) {
        for(group_b in (group_a+1):num_groups) {
          CPP_update_combinations(length_root_distances, topological_root_distances, node$partitions[[group_a]], 
                                  node$partitions[[group_b]], index_offsets, node$root_distance, node$edges_to_root)
        }
      }
      
    }
  })
  
  if(!return_lambda_function)
    return(lambda * length_root_distances + (1-lambda) * topological_root_distances)
  else {
    return(function(l) { 
      if(l<0 || l>1) stop("Pick lambda in [0,1]")
      return(l * length_root_distances + (1-l) * topological_root_distances) })
  }
}
#tree.vec <- cmpfun(tree.vec)




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
#' ## pairwise distance matrix when lambda=0
#' multi.dist(trees) 
#' ## pairwise distance matrix as a function of lambda:
#' m <- multi.dist(trees, type="function") 
#' ## evaluate at lambda=0. Equivalent to multi.dist(trees).
#' m(0)
#'
#' ## A method to visualise these distances with MDS:
#' require(ade4)
#' # find an approximate projection of the points in 2 dimensions:
#' mMDS <- dudi.pco(as.dist(m(0)), scannf=FALSE,nf=2)
#' # put the coordinates of these points into a data frame
#' mdf <- as.data.frame(cbind(mMDS$li[,1],mMDS$li[,2])) 
#' # create plot:
#' \donttest{require(ggplot2)
#' mplot <- ggplot(mdf, aes(mMDS$li[,1],mMDS$li[,2])) 
#' require(RColorBrewer)
#' mpalette <- brewer.pal(10,"Paired") # create colour palette
#' ## plot:
#' mplot + geom_point(colour=mpalette,size=5) +  
#'   xlab("") + ylab("") + theme_bw(base_family = "") }
#'
#' ## An example using data:
#' ## These woodmice phylogenies were created using the bootstrapping example in package \code{ape}
#' data(woodmiceTrees)
#' \donttest{woodmiceDists <- multi.dist(woodmiceTrees) # find topological distances
#' woodmiceMDS <- dudi.pco(as.dist(woodmiceDists), scannf=FALSE, nf=2)
#' woodmicedf <- as.data.frame(cbind(woodmiceMDS$li[,1], woodmiceMDS$li[,2]))
#' require(ggplot2)
#' woodmiceplot <- ggplot(woodmicedf, aes(woodmiceMDS$li[,1], woodmiceMDS$li[,2])) # create plot
#' woodmiceplot + geom_density2d(colour="gray80") + # contour lines
#'  geom_point(size=6, shape=1, colour="gray50") + # grey edges
#'  geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
#'  xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background}
#'
#' \donttest{require(rgl)
#' woodmiceMDS3D <- dudi.pco(as.dist(woodmiceDists), scannf=FALSE, nf=3)
#' plot3d(woodmiceMDS3D$li[,1], woodmiceMDS3D$li[,2], woodmiceMDS3D$li[,3], type="s", size=1.5, 
#'    col="navy", alpha=0.5, xlab="", ylab="", zlab="")}
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
