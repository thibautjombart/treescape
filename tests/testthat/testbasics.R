library(testthat)
library(treescape)
library(ape)

############################
# create some test objects
############################

tree_a <- rtree(100)
tree_b <- rtree(100)
trees <- rmtree(5,100)
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
  
## To add: 
# more tests for med.tree
# expect warnings and errors
