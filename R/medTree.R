
#' Geometric median tree function
#'
#' Finds the geometric median of a set of trees according to the Kendall Colijn metric.
#'
#' @export
#'
#' @author Jacob Almagro-Garcia \email{nativecoder@@gmail.com}
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @param x An object of the class multiPhylo, containing the trees for which the median tree will be computed.
#' @param groups an optional factor defining groups of trees; if provided, one median tree will be seeked for each group.
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{return.lambda.function=TRUE}.
#' @param weights A vector of weights for the trees. Defaults to a vector of 1's so that all trees are equally weighted, but can be used to encode likelihood, posterior probabilities or other characteristics.
#' @param return.lambda.function If true, a function that can be invoked with different lambda values is returned.
#'  This function returns the vector of metric values for the given lambda.
#' @param save.memory A flag that saves a lot of memory but increases the execution time (not compatible with return.lambda.function=TRUE).
#'
#' @return A list with the median metric vector, distances, indices of the tree(s) that are closest to the median tree and the value of this distance or a function that produces this list for a given value of lambda. If groups are provided, then one list is returned for each group.
#'
#'
#' @import ape
#'
#'
#' @examples
#'
#' ## generate 10 random trees, each with 6 tips
#' trees <- rmtree(10,6)
#' ## Geometric median tree:
#' mymedian <- medTree(trees)
#' mymedian$centre # the vector at the 'centre' of the trees; may not correspond to an actual tree
#' mymedian$median # the identifier(s) of the tree(s) closest to the central vector
#' mymedian$mindist # the distance of the median tree(s) from the central vector
#'
#' \dontrun{
#' ## Example with woodmice data:
#' data(woodmiceTrees)
#' woodmiceMed <- medTree(woodmiceTrees)$median[[1]]
#' ## plot the (first) geometric median tree (there are seven topologically identical median trees):
#' plot(woodmiceTrees[[woodmiceMed]],type="cladogram",edge.width=3, cex=0.8)
#'
#' ## finding the geometric median tree from a single cluster:
#' woodmiceDists <- treescape(woodmiceTrees,nf=2)
#' wmx <- woodmiceDists$pco$li[,1] # simplifying notation
#' wmy <- woodmiceDists$pco$li[,2]
#' ## isolate the trees from the largest cluster
#' wmCluster1 <- woodmiceTrees[intersect(
#'   intersect(which(wmx>(-2)),which(wmx<2)),
#'   intersect(which(wmy>(-2.5)),which(wmy<2.5))
#'   )]
#' ## find the geometric median
#' geomMedwm1 <- medTree(wmCluster1)$median[[1]]
#' plot(wmCluster1[[geomMedwm1]],type="cladogram",edge.width=3, cex=0.8)
#' # this is identical to the overall median tree:
#' treeDist(woodmiceTrees[[woodmiceMed]],wmCluster1[[geomMedwm1]],1)
#'
#' ## However, median trees from other clusters have different topologies, for example:
#' ## isolate the trees from the second largest cluster:
#' wmCluster2 <- woodmiceTrees[intersect(
#'  intersect(which(wmx>(-1)),which(wmx<8)),
#'  intersect(which(wmy>1),which(wmy<6))
#'  )]
#' ## find the geometric median
#' geomMedwm2 <- medTree(wmCluster2)$median[[1]]
#' plot(wmCluster2[[geomMedwm2]],type="cladogram",edge.width=3, cex=0.8)
#' ## This is another representative summary tree which is different from those we found above:
#' treeDist(wmCluster1[[geomMedwm1]],wmCluster2[[geomMedwm2]])
#' }
medTree <- function(x, groups=NULL, lambda=0, weights=rep(1,length(x)),
                    return.lambda.function=FALSE, save.memory=FALSE) {

    ## DEFINE MAIN FUNCTION FINDING MEDIAN TREE ##
    findMedian <- function(trees){
        num_trees <- length(trees)
        num_leaves <- length(trees[[1]]$tip.label)

        if(length(weights)!=num_trees) stop("Length of vector of weights must be the same as number of trees")

        ## Working with numbers (no functions).
        if(!return.lambda.function) {

            ## Here we speed up the computation by storing all vectors (a lot of memory for big trees).
            if(!save.memory) {

                ## Compute the metric vector for all trees.
                tree_metrics <- t(sapply(trees, function(tree) {treeVec(tree, lambda, F)}))

                ## Compute the centre metric vector by weighting the metric vector of each tree.
                centre <- (weights %*% tree_metrics)/num_trees

                ## Distances to the centre.
                distances <- apply(tree_metrics, 1, function(m){sqrt(sum((m-centre)^2))})

                ## Get the indices for the median tree(s).
                min_distance <- min(distances)
                median_trees <- which(min_distance == distances)

                return(list(centre=centre, distances=distances, mindist=min_distance, median=median_trees))
            }

            ## To save memory we recompute the vectors on the fly (way slower but we don't eat a ton of memory).
            ## We'll need a first pass to compute the centre and a second pass to compute distances.
            else {

                ## First pass: compute the centre.
                centre <- rep(0,(num_leaves*(num_leaves-1)/2) + num_leaves)
                for(i in 1:num_trees) {
                    centre <- centre + treeVec(trees[[i]], lambda, F) * weights[i]
                }
                centre <- centre/num_trees

                ## Second pass: compute the distances.
                distances <- rep(NA,num_trees)
                for(i in 1:num_trees) {
                    distances[i] <- sqrt(sum((treeVec(trees[[i]], lambda, F) - centre)^2))
                }

                ## Get the indices for the median tree(s).
                min_distance <- min(distances)
                median_trees <- which(min_distance == distances)

                return(list(centre=centre, distances=distances, mindist=min_distance, median=median_trees))
            }
        }

        ## Working with functions.
        else {

            if(save.memory)
                warning("save.memory=TRUE is incompatible with return.lambda.function=TRUE, setting save.memory=FALSE")

            ## Compute the list of metric functions for all trees.
            tree_metric_functions <- sapply(trees, function(tree) {treeVec(tree, lambda, T)})

            ## Inner function that we'll return, computes the distance matrix given lambda.
            compute_median_tree_function <- function(l) {

                ## Compute the tree metrics for the given lambda.
                tree_metrics <- t(sapply(tree_metric_functions, function(tmf){tmf(l)}))

                ## Compute the centre metric vector by weighting the metric vector of each tree.
                centre <- (weights %*% tree_metrics)/num_trees

                ## Distances to the centre.
                distances <- apply(tree_metrics, 1, function(m){sqrt(sum((m-centre)^2))})

                ## Get the indices for the median tree(s).
                min_distance <- min(distances)
                median_trees <- which(min_distance == distances)

                return(list(centre=centre, distances=distances, mindist=min_distance, median=median_trees))
            }

            return(compute_median_tree_function)
        }
    } # end findMedian


    ## APPLY FUNCTION TO TREES ##
    if(is.null(groups)){     ## no groups provided
        out <- findMedian(x)
    } else { ## groups provided
        out <- tapply(x, groups, findMedian)
    }

    ## RETURN ##
    return(out)
} ## end medTree
