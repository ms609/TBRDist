#' Tree Rearrangement Distances
#'
#'
#'
#' @param tree1,tree2 Trees of class `phylo`, or lists thereof.
#' @param checks Logical specifying whether to check that trees are properly
#' formatted and labelled.  Specify `FALSE` at your peril, as improper
#' input is likely to crash R.
#'
#' @param useApproxEstimate,useTbrEstimate,useReplugEstimate Logical specifying
#' whether to use approximate TBR distance, TBR distance or Replug distance to
#' help estimate the SPR distance.
#'
#' @name TreeRearrangementDistances
#' @export
#' @template MRS
USPRDist <- function (tree1, tree2, checks = TRUE,
                      useTbrApproxEstimate = TRUE,
                      useTbrEstimate = TRUE,
                      useReplugEstimate = TRUE) {
  treeLists <- .PrepareTrees(tree1, tree2, checks)
  uspr_dist(treeLists[[1]], treeLists[[2]], keepLabels = FALSE,
            useTbrApproxEstimate = useTbrApproxEstimate,
            useTbrEstimate = useTbrEstimate,
            useReplugEstimate = useReplugEstimate)
}

#' @rdname TreeRearrangementDistances
#' @param maf Logical specifying whether to report a maximum agreement forest
#' corresponding to the optimal score.
#' @export
ReplugDist <- function (tree1, tree2, checks = TRUE, maf = FALSE) {
  treeLists <- .PrepareTrees(tree1, tree2, checks)
  ret <- replug_dist(treeLists[[1]], treeLists[[2]], keepLabels = FALSE)
  if (maf) {
    names(ret) <- c('replug_dist', 'maf_1', 'maf_2')
    ret
  } else {
    ret[[1]]
  }
}

#' @rdname TreeRearrangementDistances
#' @param exact Logical specifying whether to calculate the exact TBR distance.
#' @param approximate Logical specifying whether to calculate the approximate
#' TBR distance.  Either this or `exact` should probably be set to `TRUE`.
#' @param countMafs Logical specifying whether to count the number of Maximum
#' Agreement Forests found.
#' @param printMafs Logical specifying whether to print Maximum Agreement
#' Forests to stdout whilst counting.
#' Use [`capture.output`]`(TBRDist(tree1, tree2, printMafs = TRUE))` to access
#' these in R.
#' @param optimize Logical specifying whether to use the default optimizations.
#' @param protectB Logical specifying whether to use the PROTECT_B optimization.
#' Overrides value inherited from `optimize`.
#' @export
TBRDist <- function (tree1, tree2, checks = TRUE,
                     exact = FALSE, approximate = !exact,
                     maf = FALSE, countMafs = FALSE, printMafs = FALSE,
                     optimize = TRUE, protectB = TRUE) {
  if (!exact && !approximate && !countMafs && !printMafs) {
    message("Nothing to do in TBRDist.")
  }
  treeLists <- .PrepareTrees(tree1, tree2, checks)

  whichRets <- c(exact, rep(approximate, 2L), countMafs,
                 rep(exact && maf, 2L))

  ret <- tbr_dist(treeLists[[1]], treeLists[[2]],
                  printMafs = printMafs, countMafs = countMafs,
                  keepLabels = FALSE,
                  optimize = optimize, protectB = protectB,
                  exact = exact, approximate = approximate)[whichRets]
  if (exact && sum(whichRets) == 1) {
    ret[[1]]
  } else {
    names(ret) <- c('tbr_exact', 'tbr_min', 'tbr_max', 'n_maf', 'maf_1', 'maf_2')[whichRets]
    ret
  }
}

#' @keywords internal
#' @importFrom ape write.tree
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
  list(write.tree(tree1), write.tree(tree2))
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
