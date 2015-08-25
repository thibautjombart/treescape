library(testthat)
library(treescape)
library(ape)

############################
# create some test objects
############################

tree_a <- rtree(100)
tree_b <- rtree(100)
n <- 10 # number of trees for multiphylo object
trees <- rmtree(n,100)
l <- runif(1) # a random value for lambda

############################
# test that evaluating at lambda immediately, or via the function, gives the same result
############################

test_that("treeVec calculated at lambda equals treeVec function evaluated at lambda", {
  expect_equal(treeVec(tree_a,l),treeVec(tree_a,return_lambda_function=TRUE)(l))
  })

test_that("treeDist calculated at lambda equals treeDist function evaluated at lambda", {
  expect_equal(treeDist(tree_a,tree_b,l),treeDist(tree_a,tree_b,return_lambda_function=TRUE)(l))
  })

test_that("multiDist calculated at lambda equals multiDist function evaluated at lambda", {
  expect_equal(multiDist(trees,l),multiDist(trees,return_lambda_function=TRUE)(l))
  })

############################
# test that functions match as they should
############################

test_that("treeDist equals Euclidean distance between corresponding vectors", {
  expect_equal(treeDist(tree_a,tree_b), sqrt(sum((treeVec(tree_a) - treeVec(tree_b))^2)))
  })

test_that("treeDist equals corresponding entry of multiDist", {
  expect_equal(treeDist(trees[[1]],trees[[2]]), multiDist(trees)[[1]])
  })
  
test_that("multiDist equals the distance matrix from treescape", {
  treedistMatrix <- treescape(trees,nf=2)$D
  expect_equal(multiDist(trees)[[n]],treedistMatrix[[n]])
  })

test_that("medTree results are consistent with treeVec", {
  geom <- medTree(trees)
  expect_equal(geom$mindist,sqrt(sum((geom$centre - treeVec(trees[[geom$median[[1]]]]))^2))) # mindist, centre and median are internally consistent, and consistent with treeVec
  expect_equal(geom$mindist,geom$distances[[geom$median[[1]]]]) # mindist equals the entry in `distances' corresponding to the (first) median tree
 })

############################
# test that save_memory versions match non-save_memory versions
############################

test_that("save_memory version of multiDist equals normal multiDist", {
  expect_equal(multiDist(trees,save_memory=TRUE), multiDist(trees))
  expect_equal(multiDist(trees,l,save_memory=TRUE), multiDist(trees,l))
  })

# NOTE: The outputs are different classes. Would like to be able to remove "as.numeric" here
test_that("save_memory version of medTree equals normal medTree", {
  expect_equal(medTree(trees,save_memory=TRUE)$centre, as.numeric(medTree(trees)$centre))
  })

############################
# test for errors and warnings
############################

test_that("error is given if lambda is outside of [0,1]", {
  expect_error(treeVec(tree_a,-1))
  expect_error(treeVec(tree_a,2))
  expect_error(treeDist(tree_a,tree_b,-1))
  expect_error(treeDist(tree_a,tree_b,2))
  expect_error(multiDist(trees,-1))
  expect_error(multiDist(trees,2))
  expect_error(medTree(trees,-1))
  expect_error(medTree(trees,2))
  })

test_that("error is given if input is not of class phylo / multiphylo", {
  expect_error(treeVec(trees))
  expect_error(treeDist(trees))
  expect_error(multiDist(tree_a))
  expect_error(medTree(tree_a))
  })

test_that("error is given if input tree is unrooted", {
  unrootedtree <- read.tree(text="(A:1,B:1,C:1);") # an unrooted tree
  expect_error(treeVec(unrootedtree))
  })

test_that("warning is given if tree edge lengths are not defined, then they are set to 1", {
  newicktree <- read.tree(text="((A,B),C);") # a tree without defined edge lengths
  expect_warning(treeVec(newicktree))
  expect_equal(treeVec(newicktree),c(1,0,0,1,1,1))
  })

test_that("error is given if trees have different tip labels", {
  tree_c <- rtree(99)
  tree_d <- tree_a
  tree_d$tip.label <- 1:100 # note that tree_a has tip labels t1, t2, ...
  expect_error(treeDist(tree_a,tree_c))
  expect_error(treeDist(tree_a,tree_d))
  })

test_that("error is given if weights vector is not of length n", {
  expect_error(medTree(trees,weights=rep(1,n+1)))
  expect_error(medTree(trees,weights=rep(1,n+1),return_lambda_function=TRUE))
  })

test_that("warning is given for the combination return_lambda_function=TRUE, save_memory=TRUE", {
  expect_warning(multiDist(trees,return_lambda_function=TRUE, save_memory=TRUE))
  expect_warning(medTree(trees,return_lambda_function=TRUE, save_memory=TRUE))
  })
