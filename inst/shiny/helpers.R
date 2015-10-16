# compute the tree vectors (as functions of lambda) only once per dataset to save on recomputation
getKCmatrixfunction <- function(x) {
  return(multiDist(x, return.lambda.function = TRUE))
}