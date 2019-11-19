#' @name Tree Rearrangement Distances
#' @export
#' @template MRS
#' @importFrom ape write.tree
USPRDist <- function (tree1, tree2, checks = TRUE) {
  treeLists <- .PrepareTrees(tree1, tree2, checks)
  uspr_dist(write.tree(treeLists[[1]]), write.tree(treeLists[[2]]), keepLabels = FALSE)
}

#' @export
#' @importFrom ape write.tree
TBRDist <- function (tree1, tree2, checks = TRUE, keepLabels = FALSE,
                     returnMaf = FALSE, printMafs = FALSE, countMafs = FALSE,
                     optimize = TRUE, protectB = TRUE,
                     exact = TRUE, approximate = !exact,
                     approxEstimate = TRUE, tbrEstimate = TRUE) {
  treeLists <- .PrepareTrees(tree1, tree2, checks)

  whichRets <- c(exact, rep(approximate, 2L), countMafs,
                 rep(exact && returnMaf, 2L))
  ret <- tbr_dist(write.tree(treeLists[[1]]), write.tree(treeLists[[2]]),
                  printMafs = printMafs, countMafs = countMafs,
                  keepLabels = FALSE,
                  optimize = optimize, protectB = protectB,
                  exact = exact, approximate = approximate,
                  approxEstimate = approxEstimate, tbrEstimate = tbrEstimate)[whichRets]
  if (exact && sum(whichRets) == 1) {
    ret[[1]]
  } else {
    names(ret) <- c('tbr_exact', 'tbr_min', 'tbr_max', 'n_maf', 'maf_1', 'maf_2')[whichRets]
    ret
  }
}

#' @keywords internal
#' @export
.PrepareTrees <- function (tree1, tree2, checks = TRUE) {
  if (checks) {

    if (class(tree1) == 'phylo') tree1 <- list(tree1)
    if (class(tree2) == 'phylo') tree2 <- list(tree2)

    if (length(tree1) != length(tree2)) {
      if (length(tree1) == 1L) {
        tree1 <- rep(tree1, length(tree2))
      } else if (length(tree2) == 1L) {
        tree2 <- rep(tree2, length(tree1))
      } else {
        stop("tree1 and tree2 must contain the same number of trees, or a single tree.")
      }
    }

    vapply(seq_along(tree1), function (i) .CatchBadPair(i, tree1[[i]], tree2[[i]]),
           integer(0))
  }
  class(tree1) <- class(tree2) <- 'multiPhylo'
  list(tree1, tree2)
}

#' @keywords internal
#' @export
.CatchBadPair <- function (i, tree1, tree2) {
  lab1 <- tree1$tip.label
  lab2 <- tree2$tip.label
  if (length(lab1) != length(lab2)) {
    stop("Problem with tree pair " , i, ": Trees must be the same size")
  }
  if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0){
    stop("Problem with tree pair " , i, ": Trees must bear identical labels")
  }
  integer(0)
}
