library(TBRDist)
library(ape)

# Regression tests for RSPRDist() against rspr 1.3.1.
# All distances were verified with the rspr command-line binary.
# Tree strings are taken from src/rspr/test_trees/ in the package source.

Tree <- function(text) suppressWarnings(read.tree(text = text))

# ── Input validation ──────────────────────────────────────────────────────────

test_that("RSPRDist rejects unrooted trees", {
  t <- unroot(rtree(6))
  expect_error(RSPRDist(t, t), "rooted")
})

test_that("RSPRDist rejects mismatched labels", {
  t1 <- rtree(6)
  t2 <- rtree(6)
  t2$tip.label <- paste0("s", seq_len(6))
  expect_error(RSPRDist(t1, t2), "identical labels")
})

# ── README examples ───────────────────────────────────────────────────────────

test_that("distance between identical trees is 0", {
  t <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  expect_identical(RSPRDist(t, t), 0L)
})

test_that("trees2.txt: exact drSPR = 4", {
  # rspr README: "total exact drSPR=4"
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));")
  expect_identical(RSPRDist(t1, t2), 4L)
})

test_that("trees2.txt: 3-approx is >= exact and <= 3x exact", {
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));")
  a <- RSPRDist(t1, t2, approx = TRUE)
  expect_gte(a, 4L)
  expect_lte(a, 12L)
})

test_that("cluster_test: exact drSPR = 3", {
  # rspr README cluster example: "total exact drSPR=3"
  t1 <- Tree("(((x,((b1,b3),b2)),y),(f,(a,c)));")
  t2 <- Tree("(((x,y),f),((a,((b1,b2),b3)),c));")
  expect_identical(RSPRDist(t1, t2), 3L)
})

# ── trees*.txt regression tests ───────────────────────────────────────────────

test_that("trees3.txt: exact drSPR = 4", {
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("((((2,1),(((11,12),4),(8,(3,(6,5))))),7),((14,(10,9)),(13,(15,16))));")
  expect_identical(RSPRDist(t1, t2), 4L)
})

test_that("trees4.txt: exact drSPR = 2", {
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("((((7,8),((2,(11,1)),(3,4))),(6,5)),((12,(10,9)),((14,13),(15,16))));")
  expect_identical(RSPRDist(t1, t2), 2L)
})

test_that("trees5.txt: exact drSPR = 1 (76 leaves)", {
  t1 <- Tree("(((((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16)))),((17,18),19)),(((((21,22),(23,24)),((25,26),(27,28))),(((29,30),(31,32)),((33,34),(35,36)))),((37,38),39))),((((((41,42),(43,44)),((45,46),(47,48))),(((49,50),(51,52)),((53,54),(55,56)))),((57,58),59)),(((((61,62),(63,64)),((65,66),(67,68))),(((69,70),(71,72)),((73,74),(75,76)))),((77,78),79))));")
  t2 <- Tree("(((((((35,36),(34,33)),((32,31),(30,29))),(((21,22),(24,23)),((25,26),(28,27)))),((38,37),39)),((18,19),((((14,13),(15,16)),((11,12),(9,10))),(((7,8),(6,5)),((3,4),(1,2)))))),((((((72,71),(70,69)),((73,74),(76,75))),(((63,64),(61,62)),((67,68),(66,65)))),(79,(77,78))),(((58,57),59),((((17,(52,51)),(50,49)),((56,55),(54,53))),(((48,47),(45,46)),((44,43),(41,42)))))));")
  expect_identical(RSPRDist(t1, t2), 1L)
})

test_that("trees6.txt: exact drSPR = 7 (39 leaves)", {
  t1 <- Tree("((((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16)))),((17,18),19)),(((((21,22),(23,24)),((25,26),(27,28))),(((29,30),(31,32)),((33,34),(35,36)))),((37,38),39)));")
  t2 <- Tree("(((((25,26),(21,22)),((31,32),(29,30))),(39,((38,((24,23),9)),37))),(((12,10),(((18,17),19),((15,16),(13,14)))),(((7,8),(5,6)),(28,(((34,33),((27,36),(35,11))),((4,3),(2,1)))))));")
  expect_identical(RSPRDist(t1, t2), 7L)
})

test_that("trees7.txt: exact drSPR = 10 (100 leaves)", {
  t1 <- Tree("((((((1,2),((3,(4,5)),6)),((7,8),(9,10))),((((((11,(12,13)),14),(15,((16,(17,(18,19))),(20,21)))),22),(23,((((24,25),26),(27,(28,29))),30))),((31,(32,33)),(((34,(35,36)),37),((38,39),((40,41),42)))))),(43,(((((44,45),46),47),((48,49),(50,((51,52),(53,54))))),((55,56),(57,58))))),((59,((60,61),(62,(63,64)))),(((65,(66,67)),((68,(69,70)),(((71,(72,(73,74))),75),(((76,77),78),(79,(80,(81,82))))))),((83,(84,((85,86),87))),((((88,(89,(90,91))),((92,93),94)),(95,(96,97))),((98,99),100))))));")
  t2 <- Tree("(((((((1,2),((3,4),6)),((7,8),(9,10))),67),((((15,((16,17),(20,21))),22),(23,(((24,26),(27,28)),30))),((31,(32,33)),((((34,36),37),((38,39),(40,42))),29)))),((43,5),(((((44,45),46),47),((48,49),(50,(((51,52),53),(18,19))))),((55,56),(57,(58,41)))))),((59,((60,61),(62,(63,64)))),(((65,66),(((68,69),((((71,((72,(73,74)),54)),70),75),(((76,77),78),(79,(80,(81,82)))))),25)),(((83,(84,((85,86),87))),((11,(12,13)),14)),((((88,(89,(90,91))),((92,93),94)),(95,(96,97))),(((98,99),100),35))))));")
  expect_identical(RSPRDist(t1, t2), 10L)
})

test_that("trees_a.txt: string labels, exact drSPR = 2", {
  t1 <- Tree("((A,B),(C,(D,E)));")
  t2 <- Tree("((A,C),(D,(B,E)));")
  expect_identical(RSPRDist(t1, t2), 2L)
})

test_that("rho_test.txt: exact drSPR = 2 (root component involved)", {
  t1 <- Tree("((132,133),(134,(135,136)));")
  t2 <- Tree("(135,(136,(134,(132,133))));")
  expect_identical(RSPRDist(t1, t2), 2L)
})

test_that("rho_test2.txt: exact drSPR = 3 (root component involved)", {
  t1 <- Tree("(((0,1),5),(7,((10,(11,12)),13)));")
  t2 <- Tree("(0,((7,((10,11),((12,5),13))),1));")
  expect_identical(RSPRDist(t1, t2), 3L)
})

# ── cluster*.txt regression tests ─────────────────────────────────────────────
# Files containing non-binary or mismatched-leaf-count pairs are omitted;
# those require rspr's -multifurcating flag, which is not exposed here.

test_that("cluster_1.txt: exact drSPR = 5", {
  t1 <- Tree("((19,(20,(21,(22,23)))),(((24,25),(26,(27,(29,((30,32),(34,35)))))),(38,(39,(40,41)))));")
  t2 <- Tree("((22,23),((26,(27,((20,19),((24,25),(21,(38,(39,(40,41)))))))),(29,((30,32),(35,34)))));")
  expect_identical(RSPRDist(t1, t2), 5L)
})

test_that("cluster_2.txt: exact drSPR = 3", {
  t1 <- Tree("((0,(1,3)),((5,(7,9)),14));")
  t2 <- Tree("(3,((7,((1,0),(5,14))),9));")
  expect_identical(RSPRDist(t1, t2), 3L)
})

test_that("cluster_3.txt: exact drSPR = 3", {
  t1 <- Tree("((((120,121),(123,129)),132),137);")
  t2 <- Tree("(((((137,132),121),129),123),120);")
  expect_identical(RSPRDist(t1, t2), 3L)
})

test_that("cluster_4.txt: exact drSPR = 12 (57 leaves, cluster decomposition)", {
  t1 <- Tree("(((17,18),(((((19,20),(38,(40,41))),(45,(46,47))),(48,(49,(50,(51,52))))),(((53,(54,55)),(56,(57,58))),((((67,((78,79),((76,77),((70,71),(72,(74,75)))))),(80,(81,(82,(83,84))))),((88,(89,(90,(91,92)))),((101,102),(94,((95,(96,97)),98))))),((60,(61,62)),(65,66)))))),((11,12),(8,9)));")
  t2 <- Tree("(55,((((56,((101,102),(94,((95,(96,97)),98)))),54),17),((88,(89,(90,(91,92)))),(((49,((67,((60,(61,62)),(65,66))),((78,79),((76,77),((70,71),(72,(74,75))))))),((53,(57,58)),((50,(51,52)),(48,((19,20),(38,(40,41))))))),((18,(((11,12),(8,9)),(80,(81,(82,(83,84)))))),(45,(46,47)))))));")
  expect_identical(RSPRDist(t1, t2), 12L)
})

test_that("cluster_6.txt: exact drSPR = 5 (multifurcating T2)", {
  t1 <- Tree("(((((0,1),2),(3,4)),(((5,6),(7,(8,(9,(10,11))))),12)),(((13,((14,(15,(16,17))),(18,(19,((20,(21,(((22,23),(24,25)),26))),(27,28)))))),(((29,30),31),32)),((33,34),35)));")
  t2 <- Tree("((3,(4,(12,((7,(8,(9,(10,11)))),(5,6))))),((32,2),((0,1),((34,33),35)),(19,((20,(27,28)),(21,(((22,23),(24,25)),26)))),(18,(14,(15,(16,17)))),13,((29,30),31)));")
  expect_identical(RSPRDist(t1, t2), 5L)
})

test_that("cluster_10.txt: exact drSPR = 1 (multifurcating T2)", {
  t1 <- Tree("(217,((218,(219,220)),221));")
  t2 <- Tree("(219,(218,217,221,220));")
  expect_identical(RSPRDist(t1, t2), 1L)
})

test_that("cluster_11.txt: exact drSPR = 7", {
  t1 <- Tree("((((((((8,(14,(26,27))),28),29),30),35),(88,(94,113))),((120,(126,127)),132)),152);")
  t2 <- Tree("(152,(127,(((30,(35,29)),((((28,14),8),27),26)),((113,126),132,120),94),88));")
  expect_identical(RSPRDist(t1, t2), 7L)
})

test_that("cluster_16.txt: exact drSPR = 2", {
  t1 <- Tree("(((0,(2,3)),(4,(6,7))),8);")
  t2 <- Tree("(((0,(3,(7,(2,4)))),6),8);")
  expect_identical(RSPRDist(t1, t2), 2L)
})

test_that("cluster_17.txt: exact drSPR = 1 (multifurcating)", {
  t1 <- Tree("(((0,1,2,3),4,(5,6,7)),(8,(9,10,((11,12,13,14),15,16,17,(18,19,20,21,22)),(23,24,25,26,27)),28,29));")
  t2 <- Tree("(((0,1,2,3,25),(5,6,7),16),((9,10,(15,17,(18,19,20,21,22,8,(11,12,13,14,4))),(23,24,26,27)),28,29));")
  expect_identical(RSPRDist(t1, t2), 1L)
})

# ── MAF output ────────────────────────────────────────────────────────────────

test_that("maf = TRUE returns correct structure and distance", {
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));")
  result <- RSPRDist(t1, t2, maf = TRUE)
  expect_named(result, c("exact", "maf_1", "maf_2"))
  expect_identical(result$exact, 4L)
  expect_gt(nchar(result$maf_1), 0L)
  expect_gt(nchar(result$maf_2), 0L)
})

test_that("MAF has exactly drSPR + 1 components", {
  # For trees_a.txt: drSPR = 2, so MAF should have 3 components
  t1 <- Tree("((A,B),(C,(D,E)));")
  t2 <- Tree("((A,C),(D,(B,E)));")
  result <- RSPRDist(t1, t2, maf = TRUE)
  n_components <- length(strsplit(trimws(result$maf_1), "\\s+")[[1]])
  expect_identical(n_components, result$exact + 1L)
})

# ── allPairs and dist object ──────────────────────────────────────────────────

test_that("allPairs = TRUE returns a dist object", {
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));")
  d <- RSPRDist(list(t1, t2))
  expect_s3_class(d, "dist")
  expect_equal(as.vector(d), 4L)
})

test_that("allPairs across two lists returns a matrix", {
  t1 <- Tree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))));")
  t2 <- Tree("(((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)));")
  m <- RSPRDist(list(t1, t2), list(t2, t1), allPairs = TRUE)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 2L))
})
