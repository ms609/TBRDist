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
  treeR <- ape::read.tree(text="(t1, (((t5, t7), (t9, (t3, t2))), (t4, ((t6, t8), t10))));")
  list1 <- list(one = tree1, oneAgain = tree1, two = tree2, three = treeR)
  list2 <- list(tree1, tree2)
  expect_equivalent(2L, USPRDist(tree1, tree2))
  expect_equivalent(c(0, 2L), USPRDist(list(tree1, tree2), tree1))
  expect_equivalent(c(0, 2L), USPRDist(tree1, list(tree1, tree2)))

  goodRet <- structure(c(0, 2, 5, 2, 5, 4),
                       Size = 4L,
                       Labels = names(list1),
                       Diag = FALSE,
                       Upper = FALSE,
                       class = 'dist')
  expect_equal(goodRet, USPRDist(list1))

  expect_equal(goodRet, ReplugDist(list1))
  first <- ReplugDist(list1, list1[[1]], maf = TRUE)
  each <- ReplugDist(list1, maf = TRUE)
  expect_equivalent(first[[1]], as.matrix(each[[1]])[, 1])
  expect_equivalent(first[[2]][-1], each[[2]][-1, 1])
  expect_equivalent(first[[3]][-1], each[[3]][-1, 1])

  TBRDist(list1, maf = TRUE)


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

