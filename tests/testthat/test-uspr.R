context("uspr.R")

test_that("SPR distances are calculated correctly", {
  tree1 <- TreeTools::BalancedTree(10)
  tree2 <- TreeTools::PectinateTree(10)
  expect_equal(phangorn::SPR.dist(tree1, tree2),
               USPRDist(tree1, tree2))

})
