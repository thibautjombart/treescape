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

test_that("tree.vec calculated at lambda equals tree.vec function evaluated at lambda", {
  expect_equal(tree.vec(tree_a,l),tree.vec(tree_a,return_lambda_function=TRUE)(l))
  })

test_that("tree.dist calculated at lambda equals tree.dist function evaluated at lambda", {
  expect_equal(tree.dist(tree_a,tree_b,l),tree.dist(tree_a,tree_b,return_lambda_function=TRUE)(l))
  })  

test_that("multi.dist calculated at lambda equals multi.dist function evaluated at lambda", {
  expect_equal(multi.dist(trees,l),multi.dist(trees,return_lambda_function=TRUE)(l))
  })    
  
############################
# test that functions match as they should
############################
  
test_that("tree.dist equals Euclidean distance between corresponding vectors", {
  expect_equal(tree.dist(tree_a,tree_b), sqrt(sum((tree.vec(tree_a) - tree.vec(tree_b))^2)))
  })

test_that("tree.dist equals corresponding entry of multi.dist", {
  expect_equal(tree.dist(trees[[1]],trees[[2]]), multi.dist(trees)[[1]])
  })  
  
test_that("med.tree results are consistent with tree.vec", {
  geom <- med.tree(trees)
  expect_equal(geom$mindist,sqrt(sum((geom$centre - tree.vec(trees[[geom$median[[1]]]]))^2))) # mindist, centre and median are internally consistent, and consistent with tree.vec
  expect_equal(geom$mindist,geom$distances[[geom$median[[1]]]]) # mindist equals the entry in `distances' corresponding to the (first) median tree
 })

############################
# test that save_memory versions match non-save_memory versions
############################

test_that("save_memory version of multi.dist equals normal multi.dist", {
  expect_equal(multi.dist(trees,save_memory=TRUE), multi.dist(trees))
  expect_equal(multi.dist(trees,l,save_memory=TRUE), multi.dist(trees,l))
  })  

# NOTE: The outputs are different classes. Would like to be able to remove "as.numeric" here
test_that("save_memory version of med.tree equals normal med.tree", {
  expect_equal(med.tree(trees,save_memory=TRUE)$centre, as.numeric(med.tree(trees)$centre))
  })  
  
############################
# test for errors and warnings
############################

test_that("error is given if lambda is outside of [0,1]", {
  expect_error(tree.vec(tree_a,-1))
  expect_error(tree.vec(tree_a,2))
  expect_error(tree.dist(tree_a,tree_b,-1))
  expect_error(tree.dist(tree_a,tree_b,2))
  expect_error(multi.dist(trees,-1))
  expect_error(multi.dist(trees,2))
  expect_error(med.tree(trees,-1))
  expect_error(med.tree(trees,2))
  })

test_that("warning is given for the combination return_lambda_function=TRUE, save_memory=TRUE", {
  expect_warning(multi.dist(trees,return_lambda_function=TRUE, save_memory=TRUE))
  expect_warning(med.tree(trees,return_lambda_function=TRUE, save_memory=TRUE))
  })
  
test_that("error is given if input is not of class phylo / multiphylo", {
  expect_error(tree.vec(trees))
  expect_error(tree.dist(trees))
  expect_error(multi.dist(tree_a))
  expect_error(med.tree(tree_a))
  })

test_that("warning is given if tree edge lengths are not defined, then they are set to 1", {
  newicktree <- read.tree(text="((A,B),C);") # a tree without defined edge lengths
  expect_warning(tree.vec(newicktree))
  expect_equal(tree.vec(newicktree),c(1,0,0,1,1,1))
  })  
  
test_that("error is given if trees have different tip labels", {
  tree_c <- rtree(99)
  tree_d <- tree_a
  tree_d$tip.label <- 1:100 # note that tree_a has tip labels t1, t2, ...
  expect_error(tree.dist(tree_a,tree_c))
  expect_error(tree.dist(tree_a,tree_d))
  })
  
test_that("error is given if weights vector is not of length n", {
  expect_error(med.tree(trees,weights=rep(1,n+1)))
  expect_error(med.tree(trees,weights=rep(1,n+1),return_lambda_function=TRUE))
  })