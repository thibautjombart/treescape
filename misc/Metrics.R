require(ape)
require(phangorn)
require(compiler)
require(fastmatch)
require(combinat)

# function to make mrca matrix of a tree, where entry (i,j) gives the mrca of tips i and j
my.mrca <- function(tree,k)
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
my.mrca <- cmpfun(my.mrca) # compile

# function to create a vector of the pendant edges of the tree
edge.pendant <- function(tree,k) {tree$edge[fmatch(1:k, tree$edge[,2]),] }
edge.pendant <- cmpfun(edge.pendant)

# function to take two objects of class phylo and return their topological distance
CK.gdist <- function(tr1,tr2) {
  k <- length(tr1$tip.label)
  # checks and warnings
  if (k != length(tr2$tip.label)) {
    stop("trees have different numbers of tips")
    }
  
  if (setequal(tr1$tip.label,tr2$tip.label) == FALSE) {
    stop("trees have different tip label sets")
    } 
  
  # set all edge lengths to 1:
  tr1$edge.length <- rep(1,(2*k-2));
  tr2$edge.length <- tr1$edge.length;
  # create vector of distances from root to tips:
  DN1 <- dist.nodes(tr1)[k+1,];
  DN2 <- dist.nodes(tr2)[k+1,];
  # create mrca matrix for each tree:
  M1 <- my.mrca(tr1,k);
  M2 <- my.mrca(tr2,k);
  # find the permutation which maps tr1 tip labels to tr2 tip labels
  labelmatch <- fmatch(tr1$tip.label, tr2$tip.label);
  
  # calculate the gamma vector of root to mrca distances:
  gamma <- apply(combn2(1:k), 1, function(x) DN1[M1[[x[1],x[2]]]] - DN2[M2[[labelmatch[x[1]],labelmatch[x[2]]]]]) 
  # Euclidean norm:
  d <- sqrt(sum(gamma^2));
  return(d)
}
# compile:
CK.gdist <- cmpfun(CK.gdist)

# function to take two objects of class phylo and return their weighted distance
CK.wdist <- function(tr1,tr2,p=1) {
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
  
  tr1$edge.length <- tr1$edge.length + p;
  tr2$edge.length <- tr2$edge.length + p;
  DN1 <- dist.nodes(tr1);
  DN2 <- dist.nodes(tr2);
  M1 <- my.mrca(tr1,k);
  M2 <- my.mrca(tr2,k);
  labelmatch <- fmatch(tr1$tip.label, tr2$tip.label);
  pairs <-combn2(1:k);
  
  delta <- apply(pairs, 1, function(x) DN1[k+1,x[1]]-DN1[k+1,x[2]] - DN2[k+1,labelmatch[x[1]]] + DN2[k+1,labelmatch[x[2]]])
  
  gamma <- apply(pairs, 1, function(x) DN1[k+1,M1[[x[1],x[2]]]] - DN2[k+1,M2[[labelmatch[x[1]],labelmatch[x[2]]]]]) 
  
  d1 <-sqrt(sum(delta^2));
  d2 <- sqrt(sum(gamma^2));
  
  # add difference in minimum pendant branch length
  ep1 <- edge.pendant(tr1,k);
  ep2 <- edge.pendant(tr2,k);
  pen.length1 <- apply(ep1, 1, function(x) DN1[x[1],x[2]])
  pen.length2 <- apply(ep2, 1, function(x) DN2[x[1],x[2]])
  near1 <- min(pen.length1);
  near2 <- min(pen.length2);
  
  d1+d2+(abs(near1-near2)/k)
}
CK.wdist <- cmpfun(CK.wdist)

# function to take an object of class multiPhylo and return topological distance matrix 
CK.gdistm <- function(trees) {
  # checks and warnings
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
  
  # set all edge lengths to 1:
  for (i in 1:l) {
    trees[[i]]$edge.length <- rep(1,(2*k-2));
  }
  DN <- sapply(1:l, function(x) dist.nodes(trees[[x]])[k+1,]); # vector of vectors
  M <- lapply(1:l, function(x) my.mrca(trees[[x]],k)); # list of matrices
  labelmatch <- lapply(1:l, function (y) fmatch(trees[[1]]$tip.label,trees[[y]]$tip.label)); # list of permutations
  pairs <-combn2(1:k);
  
  gammas <- sapply(1:l, function(y) apply(pairs, 1, function(x) DN[M[[y]][[labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]],y])); 
  
  x <- k*(k-1)/2;
  
  # create lower triangular matrix
  # first, where each entry is a vector of differences
  gammadistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
              if (a<b) sapply(1:x, function(c) (gammas[[c,a]]-gammas[[c,b]])^2))) 
  # now find Euclidean norm of each entry        
  gammadistmat <- sapply(1:l, function(a) 
                      sapply(1:l, function(b) sqrt(sum(gammadistmat[[b,a]]))))
  
  return(gammadistmat)
}
CK.gdistm <- cmpfun(CK.gdistm)


# function to take an object of class multiPhylo and return weighted distance matrix 
CK.wdistm <- function(trees,p=1) {
  # checks and warnings
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
    if (is.null(trees[[i]]$edge.length)) {
      stop("edge lengths not defined")
    }
    trees[[i]]$edge.length <- trees[[i]]$edge.length + p;
  }
  
  DN <- lapply(1:l, function(x) dist.nodes(trees[[x]]));
  M <- lapply(1:l, function(x) my.mrca(trees[[x]],k));
  labelmatch <- lapply(1:l, function (y) 
                  fmatch(trees[[1]]$tip.label,trees[[y]]$tip.label));
  pairs <-combn2(1:k);
  
  deltas <- sapply(1:l, function(y) apply(pairs, 1, function(x) 
    DN[[y]][k+1,labelmatch[[y]][x[1]]] - DN[[y]][k+1,labelmatch[[y]][x[2]]]))
  gammas <- sapply(1:l, function(y) apply(pairs, 1, function(x) 
    DN[[y]][k+1,M[[y]][labelmatch[[y]][x[1]],labelmatch[[y]][x[2]]]])); 
  
  x <- k*(k-1)/2;
  
  deltadistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
    if (a<b) sapply(1:x, function(c) (deltas[[c,a]]-deltas[[c,b]])^2)))
  deltadistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
                        sqrt(sum(deltadistmat[[b,a]]))))
  
  gammadistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
    if (a<b) sapply(1:x, function(c) (gammas[[c,a]]-gammas[[c,b]])^2)))
  gammadistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
                        sqrt(sum(gammadistmat[[b,a]]))))
  
  E <- lapply(1:l, function(x) edge.pendant(trees[[x]],k))
  P <- sapply(1:l, function(x) apply(E[[x]], 1, function(y) DN[[x]][y[1],y[2]]) )
  minP <- sapply(1:l, function(x) min(P[,x]))
  sigmadistmat <- sapply(1:l, function(a) sapply(1:l, function(b) 
    if (a<b) {abs(minP[[a]]-minP[[b]])/k} else {0}))
  
  distmat <- deltadistmat + gammadistmat + sigmadistmat
  return(distmat)
}
CK.wdistm <- cmpfun(CK.wdistm)