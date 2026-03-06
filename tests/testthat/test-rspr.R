library(TBRDist)
library(ape)

# Reference distances verified against rspr 1.3.1 command-line binary.

# ── Helper ────────────────────────────────────────────────────────────────────

newick_pair <- function(t1, t2) {
  list(read.tree(text = t1), read.tree(text = t2))
}

# ── Input validation ──────────────────────────────────────────────────────────

test_that("RSPRDist rejects unrooted trees", {
  t <- unroot(rtree(6))
  expect_error(RSPRDist(t, t), "rooted")
})

test_that("RSPRDist rejects mismatched labels", {
  t1 <- rtree(6)
  t2 <- rtree(6)
  t2$tip.label <- paste0("s", seq_len(6))  # ensure different labels
  expect_error(RSPRDist(t1, t2), "identical labels")
})

# ── Distances ─────────────────────────────────────────────────────────────────

test_that("distance between identical trees is 0", {
  tr <- newick_pair(
    "((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));",
    "((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));"
  )
  expect_identical(RSPRDist(tr[[1]], tr[[2]]), 0L)
})

test_that("rspr exact: trees2.txt pair gives distance 4", {
  # Reference: rspr README / example output: "total exact drSPR=4"
  tr <- newick_pair(
    "((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));",
    "(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));"
  )
  expect_identical(RSPRDist(tr[[1]], tr[[2]]), 4L)
})

test_that("rspr approx: trees2.txt pair gives at most 3x exact", {
  tr <- newick_pair(
    "((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));",
    "(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));"
  )
  # Exact is 4; 3-approx must be >= exact and <= 3 * exact
  a <- RSPRDist(tr[[1]], tr[[2]], approx = TRUE)
  expect_gte(a, 4L)
  expect_lte(a, 12L)
})

test_that("rspr exact: README cluster example gives distance 3", {
  # Reference: rspr README: "total exact drSPR=3"
  tr <- newick_pair(
    "(((x,((b1,b3),b2)),y),(f,(a,c)));",
    "(((x,y),f),((a,((b1,b2),b3)),c));"
  )
  expect_identical(RSPRDist(tr[[1]], tr[[2]]), 3L)
})

# ── MAF output ────────────────────────────────────────────────────────────────

test_that("maf = TRUE returns a list with correct structure", {
  tr <- newick_pair(
    "((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));",
    "(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));"
  )
  result <- RSPRDist(tr[[1]], tr[[2]], maf = TRUE)
  expect_named(result, c("exact", "maf_1", "maf_2"))
  expect_identical(result$exact, 4L)
  # MAFs are non-empty strings
  expect_gt(nchar(result$maf_1), 0L)
  expect_gt(nchar(result$maf_2), 0L)
})

# ── allPairs and dist object ──────────────────────────────────────────────────

test_that("allPairs returns a dist object for a list of trees", {
  t1 <- read.tree(text = "((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- read.tree(text = "(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));")
  trees <- list(t1, t2)
  d <- RSPRDist(trees)
  expect_s3_class(d, "dist")
  expect_equal(as.vector(d), 4L)
})
