context("uspr.R")
library('TreeTools')

test_that("Bad tree input is handled correctly", {
  expect_error(USPRDist(PectinateTree(1:8), PectinateTree(2:9)))
  expect_error(USPRDist(PectinateTree(1:8), PectinateTree(1:9)))
  expect_error(USPRDist(PectinateTree(1:9), PectinateTree(1:8)))
  expect_error(USPRDist(list(PectinateTree(1:8), BalancedTree(1:8)),
                        list(PectinateTree(1:8), BalancedTree(1:8), BalancedTree(1:8))))
})

test_that("SPR distances are calculated correctly", {
  tree1 <- BalancedTree(10)
  tree2 <- PectinateTree(10)
  expect_equivalent(2L, USPRDist(tree1, tree2))
  expect_equivalent(c(0, 2L), USPRDist(list(tree1, tree2), tree1))
  expect_equivalent(c(0, 2L), USPRDist(tree1, list(tree1, tree2)))


  Test <- function (tree1, tree2) {
    td <- TBRDist(tree1, tree2, exact = TRUE, approximate = TRUE)
    expect_true(USPRDist(tree1, tree2) >= TBRDist(tree1, tree2, exact = TRUE))
    expect_true(td$tbr_exact >= td$tbr_min)
    expect_true(td$tbr_exact <= td$tbr_max)

    td4 <- TBRDist(list(tree1, tree1, tree2, tree2),
                   list(tree1, tree2, tree1, tree2), exact = TRUE)
    expect_equal(0L, td4[[1]])
    expect_equal(td4[[2]], td4[[3]])
    expect_equal(0L, td4[[4]])

    sd4 <- USPRDist(list(tree1, tree1, tree2, tree2),
                    list(tree1, tree2, tree1, tree2))
    expect_equal(0L, sd4[[1]])
    expect_equal(sd4[[2]], sd4[[3]])
    expect_equal(0L, sd4[[4]])

    rd4 <- ReplugDist(list(tree1, tree1, tree2, tree2),
                      list(tree1, tree2, tree1, tree2))
    expect_equal(0L, rd4[[1]])
    expect_equal(rd4[[2]], rd4[[3]])
    expect_equal(0L, rd4[[4]])

  }

  Test(tree1, tree2)
  Test(PectinateTree(13), BalancedTree(13))

  expect_equal(invisible(),
               TBRDist(tree1, tree2, exact = FALSE, approximate = FALSE))
})
