context("uspr.R")
library('TreeTools')

test_that("Bad tree input is handled correctly", {
  expect_error(USPRDist(PectinateTree(1:8), PectinateTree(2:9)))
  expect_error(USPRDist(PectinateTree(1:8), PectinateTree(1:9)))
  expect_error(USPRDist(PectinateTree(1:9), PectinateTree(1:8)))

})

test_that("SPR distances are calculated correctly", {
  tree1 <- TreeTools::BalancedTree(10)
  tree2 <- TreeTools::PectinateTree(10)
  expect_equivalent(phangorn::SPR.dist(tree1, tree2),
                    USPRDist(tree1, tree2))

  expect_true(USPRDist(tree1, tree2) >= TBRDist(tree1, tree2))
})
